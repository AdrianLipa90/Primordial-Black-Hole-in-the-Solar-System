#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import io
import json
import math
import sys
from typing import Any

import numpy as np
import requests
from astropy.io import fits

TARGET = {
    "ra_deg": 107.576646,
    "dec_deg": -25.907579,
    "gal_l_deg": 238.1709218507047,
    "gal_b_deg": -7.665782897020605,
}
TILE_ORDER = 3
TILE_NPIX = 488

SOURCES = {
    "planck_100": {
        "url": "https://alasky.cds.unistra.fr/PLANCK/R3/HFI_SkyMap_100_R3/Norder3/Dir0/Npix488.fits",
        "kind": "CMB_FREQUENCY",
    },
    "planck_143": {
        "url": "https://alasky.cds.unistra.fr/PLANCK/R3/HFI_SkyMap_143_R3/Norder3/Dir0/Npix488.fits",
        "kind": "CMB_FREQUENCY",
    },
    "planck_217": {
        "url": "https://alasky.cds.unistra.fr/PLANCK/R3/HFI_SkyMap_217_R3/Norder3/Dir0/Npix488.fits",
        "kind": "CMB_FREQUENCY",
    },
    "wmap_V": {
        "url": "https://alasky.cds.unistra.fr/WMAP9yr/WMAPV9yr/Norder3/Dir0/Npix488.fits",
        "kind": "CMB_FREQUENCY",
    },
    "wmap_W": {
        "url": "https://alasky.cds.unistra.fr/WMAP9yr/WMAPW9yr/Norder3/Dir0/Npix488.fits",
        "kind": "CMB_FREQUENCY",
    },
}


def spread_bits(value: int) -> int:
    out = 0
    bit = 0
    while value:
        if value & 1:
            out |= 1 << (2 * bit)
        value >>= 1
        bit += 1
    return out


def compact_bits(value: int) -> int:
    out = 0
    bit = 0
    src = 0
    while value >> src:
        if value & (1 << src):
            out |= 1 << bit
        bit += 1
        src += 2
    return out


def ang2pix_nest(nside: int, lon_deg: float, lat_deg: float) -> int:
    if nside <= 0 or nside & (nside - 1):
        raise ValueError("nside must be a positive power of two")

    phi = math.radians(lon_deg) % (2.0 * math.pi)
    z = math.sin(math.radians(lat_deg))
    za = abs(z)
    tt = phi / (0.5 * math.pi)

    if za <= 2.0 / 3.0:
        temp1 = nside * (0.5 + tt)
        temp2 = nside * (0.75 * z)
        jp = math.floor(temp1 - temp2)
        jm = math.floor(temp1 + temp2)
        ifp = jp // nside
        ifm = jm // nside
        if ifp == ifm:
            face = ifp | 4
        elif ifp < ifm:
            face = ifp
        else:
            face = ifm + 8
        ix = jm % nside
        iy = nside - (jp % nside) - 1
    else:
        ntt = math.floor(tt)
        tp = tt - ntt
        tmp = nside * math.sqrt(3.0 * (1.0 - za))
        jp = min(math.floor(tp * tmp), nside - 1)
        jm = min(math.floor((1.0 - tp) * tmp), nside - 1)
        if z >= 0:
            face = ntt
            ix = nside - jm - 1
            iy = nside - jp - 1
        else:
            face = ntt + 8
            ix = jp
            iy = jm

    return int(face * nside * nside + spread_bits(ix) + (spread_bits(iy) << 1))


def local_tile_xy(tile_width: int) -> tuple[int, int, int]:
    bits = int(round(math.log2(tile_width)))
    if 2**bits != tile_width:
        raise ValueError("HiPS tile width must be a power of two")
    final_order = TILE_ORDER + bits
    nside = 2**final_order
    pix = ang2pix_nest(nside, TARGET["gal_l_deg"], TARGET["gal_b_deg"])
    parent = pix >> (2 * bits)
    if parent != TILE_NPIX:
        raise RuntimeError(
            f"target child pixel parent {parent} != preregistered tile {TILE_NPIX}"
        )
    local = pix & ((1 << (2 * bits)) - 1)
    x = compact_bits(local)
    y = compact_bits(local >> 1)
    return x, y, final_order


def finite_image(image: np.ndarray) -> np.ndarray:
    image = np.asarray(image, dtype=float)
    finite = np.isfinite(image)
    if not finite.any():
        raise RuntimeError("image contains no finite pixels")
    med = float(np.nanmedian(image[finite]))
    return np.where(finite, image, med)


def phase_markers(image: np.ndarray) -> dict[str, float | list[int] | None]:
    a = finite_image(image)
    ny, nx = a.shape
    if min(nx, ny) < 16:
        return {"error": "PATCH_TOO_SMALL"}

    raw_mean = float(np.mean(a))
    raw_std = float(np.std(a))
    centered = a - float(np.median(a))
    window = np.outer(np.hanning(ny), np.hanning(nx))
    fourier = np.fft.fftshift(np.fft.fft2(centered * window))
    amp = np.abs(fourier)
    power = amp * amp
    phase = np.angle(fourier)

    yy, xx = np.indices(a.shape)
    cy, cx = (ny - 1) / 2.0, (nx - 1) / 2.0
    kx = (xx - cx) / nx
    ky = (yy - cy) / ny
    kr = np.hypot(kx, ky)
    theta = np.arctan2(ky, kx)

    band = (kr > 0.03) & (kr < 0.22)
    aw = amp[band]
    pw = power[band]
    if not band.any() or float(np.sum(aw)) == 0.0:
        return {"error": "EMPTY_PHASE_BAND"}

    coherence = float(
        abs(np.sum(aw * np.exp(1j * phase[band]))) / np.sum(aw)
    )

    hist, _ = np.histogram(
        phase[band], bins=24, range=(-np.pi, np.pi), weights=aw
    )
    p = hist / np.sum(hist)
    p = p[p > 0]
    entropy = float(-np.sum(p * np.log(p)) / np.log(24))

    out: dict[str, Any] = {
        "shape": [int(ny), int(nx)],
        "mean": raw_mean,
        "std": raw_std,
        "phase_coherence": coherence,
        "phase_entropy": entropy,
    }
    for n in (1, 2, 3, 6):
        out[f"A{n}"] = float(
            abs(np.sum(pw * np.exp(1j * n * theta[band]))) / np.sum(pw)
        )
    return out


def target_patch(data: np.ndarray, x: int, y: int) -> tuple[np.ndarray, dict[str, Any]]:
    ny, nx = data.shape
    # HiPS image y-axis orientation can differ by writer/viewer convention.
    candidates = {
        "direct": (x, y),
        "y_flipped": (x, ny - 1 - y),
    }
    half = max(8, min(nx, ny) // 8)
    patches: dict[str, Any] = {}
    for label, (px, py) in candidates.items():
        x0, x1 = max(0, px - half), min(nx, px + half + 1)
        y0, y1 = max(0, py - half), min(ny, py + half + 1)
        patch = data[y0:y1, x0:x1]
        patches[label] = {
            "xy": [int(px), int(py)],
            "bounds": [int(x0), int(x1), int(y0), int(y1)],
            "markers": phase_markers(patch),
            "center_value": float(data[py, px]),
        }
    # Return direct for pairwise matching, while preserving both orientation
    # receipts so no image-axis assumption is hidden.
    px, py = candidates["direct"]
    x0, x1 = max(0, px - half), min(nx, px + half + 1)
    y0, y1 = max(0, py - half), min(ny, py + half + 1)
    return data[y0:y1, x0:x1], patches


def downsample_to(image: np.ndarray, target: int) -> np.ndarray:
    image = finite_image(image)
    ny, nx = image.shape
    if ny != nx:
        raise ValueError("only square HiPS tiles are supported")
    if nx == target:
        return image
    if nx % target != 0:
        raise ValueError(f"cannot block-average {nx} to {target}")
    factor = nx // target
    return image.reshape(target, factor, target, factor).mean(axis=(1, 3))


def corr(a: np.ndarray, b: np.ndarray) -> float | None:
    size = min(a.shape[0], b.shape[0])
    a2 = downsample_to(a, size)
    b2 = downsample_to(b, size)
    x = (a2 - np.mean(a2)).ravel()
    y = (b2 - np.mean(b2)).ravel()
    denom = float(np.linalg.norm(x) * np.linalg.norm(y))
    return float(np.dot(x, y) / denom) if denom > 0 else None


def main() -> int:
    receipt: dict[str, Any] = {
        "schema": "PBH_CMB_HIPS_PHASE_RECEIPT_V0_1",
        "target": TARGET,
        "tile": {"order": TILE_ORDER, "npix": TILE_NPIX},
        "maps": {},
        "pairwise_tile_correlations": {},
        "status": "PARTIAL",
        "interpretation": "DESCRIPTIVE_ONLY_NOT_PBH_DETECTION",
        "authority": "CANDIDATE_ONLY",
        "canon_allowed": False,
    }
    images: dict[str, np.ndarray] = {}

    for name, source in SOURCES.items():
        entry: dict[str, Any] = {
            "source_url": source["url"],
            "kind": source["kind"],
        }
        try:
            response = requests.get(
                source["url"],
                timeout=120,
                headers={"User-Agent": "GREMLIN-FERRET-PBH-CMB/0.1 read-only"},
            )
            entry["http_status"] = response.status_code
            entry["bytes"] = len(response.content)
            entry["sha256"] = hashlib.sha256(response.content).hexdigest()
            entry["content_type"] = response.headers.get("content-type")
            if response.status_code != 200:
                raise RuntimeError(f"HTTP {response.status_code}")
            if not response.content.startswith(b"SIMPLE"):
                raise RuntimeError("response is not a FITS primary HDU")

            with fits.open(io.BytesIO(response.content), memmap=False) as hdul:
                hdu = hdul[0]
                data = np.asarray(hdu.data, dtype=float)
                while data.ndim > 2:
                    data = data[0]
                if data.ndim != 2:
                    raise RuntimeError(f"expected 2D FITS tile, got {data.shape}")
                data = finite_image(data)
                header = hdu.header

            tile_width = int(data.shape[1])
            if data.shape[0] != tile_width:
                raise RuntimeError(f"non-square HiPS tile: {data.shape}")
            x, y, final_order = local_tile_xy(tile_width)
            patch, target_receipt = target_patch(data, x, y)

            entry.update(
                {
                    "shape": [int(data.shape[0]), int(data.shape[1])],
                    "tile_width": tile_width,
                    "final_healpix_order": final_order,
                    "target_local_xy_nested": [x, y],
                    "header": {
                        key: header.get(key)
                        for key in (
                            "BUNIT",
                            "BSCALE",
                            "BZERO",
                            "OBJECT",
                            "EXTNAME",
                            "NAXIS1",
                            "NAXIS2",
                        )
                        if key in header
                    },
                    "global_markers": phase_markers(data),
                    "target_orientation_receipts": target_receipt,
                }
            )
            images[name] = data
            receipt["maps"][name] = entry
        except Exception as exc:
            entry["error"] = f"{type(exc).__name__}: {exc}"
            receipt["maps"][name] = entry

    names = sorted(images)
    for i, left in enumerate(names):
        for right in names[i + 1 :]:
            try:
                receipt["pairwise_tile_correlations"][f"{left}__{right}"] = corr(
                    images[left], images[right]
                )
            except Exception as exc:
                receipt["pairwise_tile_correlations"][f"{left}__{right}"] = {
                    "error": f"{type(exc).__name__}: {exc}"
                }

    ok = [k for k, v in receipt["maps"].items() if "error" not in v]
    receipt["successful_maps"] = ok
    receipt["status"] = "ANALYZED" if len(ok) == len(SOURCES) else "PARTIAL"

    payload = json.dumps(receipt, sort_keys=True, separators=(",", ":"))
    print("PBH_CMB_RECEIPT=" + payload)
    return 0 if ok else 2


if __name__ == "__main__":
    sys.exit(main())
