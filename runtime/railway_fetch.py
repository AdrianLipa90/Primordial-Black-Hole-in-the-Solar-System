#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import io
import json
import math
import sys
from typing import Any
from urllib.parse import urlencode

import numpy as np
import requests
from astropy.io import fits

TARGET = {
    "ra_deg": 107.576646,
    "dec_deg": -25.907579,
    "gal_l_deg": 238.1709218507047,
    "gal_b_deg": -7.665782897020605,
}

WIDTH = 256
HEIGHT = 256
FOV_DEG = 8.0
PIX_ARCMIN = FOV_DEG * 60.0 / WIDTH
COMMON_FWHM_ARCMIN = 40.0

SOURCES = {
    "planck_smica": {
        "hips": "CDS/P/PLANCK/R3/CMB",
        "kind": "CMB_COMPONENT",
        "native_fwhm_arcmin": 5.0,
    },
    "planck_100": {
        "hips": "CDS/P/PLANCK/R3/HFI100",
        "kind": "CMB_FREQUENCY",
        "native_fwhm_arcmin": 9.66,
    },
    "planck_143": {
        "hips": "CDS/P/PLANCK/R3/HFI143",
        "kind": "CMB_FREQUENCY",
        "native_fwhm_arcmin": 7.27,
    },
    "planck_217": {
        "hips": "CDS/P/PLANCK/R3/HFI217",
        "kind": "CMB_FREQUENCY",
        "native_fwhm_arcmin": 5.01,
    },
    "wmap_V": {
        "hips": "CDS/P/WMAP/V/9yr",
        "kind": "CMB_FREQUENCY",
        "native_fwhm_arcmin": 21.0,
    },
    "wmap_W": {
        "hips": "CDS/P/WMAP/W/9yr",
        "kind": "CMB_FREQUENCY",
        "native_fwhm_arcmin": 13.2,
    },
}

HIPS2FITS = "https://alasky.cds.unistra.fr/hips-image-services/hips2fits"


def finite_image(image: np.ndarray) -> np.ndarray:
    image = np.asarray(image, dtype=float)
    finite = np.isfinite(image)
    if not finite.any():
        raise RuntimeError("image contains no finite pixels")
    med = float(np.nanmedian(image[finite]))
    return np.where(finite, image, med)


def phase_markers(image: np.ndarray) -> dict[str, Any]:
    a = finite_image(image)
    ny, nx = a.shape
    if min(nx, ny) < 32:
        return {"error": "PATCH_TOO_SMALL"}

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

    band = (kr > 0.02) & (kr < 0.20)
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
        "mean": float(np.mean(a)),
        "std": float(np.std(a)),
        "phase_coherence": coherence,
        "phase_entropy": entropy,
    }
    for n in (1, 2, 3, 6):
        out[f"A{n}"] = float(
            abs(np.sum(pw * np.exp(1j * n * theta[band]))) / np.sum(pw)
        )
    return out


GREMLIN_SOURCE = {
    "repository": "AdrianLipa90/GREMLIN",
    "commit": "9fdc42e507d5fe54c13dda975c297173ba0e4eb6",
    "artifact": "gremlin_mcp/holographic_orbital_memory.py",
    "algebra": [
        "36D phase cell",
        "phase_rms_distance",
        "moire_overlap",
        "project_hologram",
    ],
    "authority": "CANDIDATE_ONLY",
}


def phase36_sector_cell(image: np.ndarray) -> dict[str, Any]:
    """Measurement-derived T^36 cell: 6 radial x 6 angular Fourier sectors.

    No semantic meanings are assigned to the 36 coordinates. Each coordinate is
    the amplitude-weighted circular mean Fourier phase of one fixed sector.
    """
    a = finite_image(image)
    ny, nx = a.shape
    centered = a - float(np.median(a))
    window = np.outer(np.hanning(ny), np.hanning(nx))
    fourier = np.fft.fftshift(np.fft.fft2(centered * window))
    amp = np.abs(fourier)
    phase = np.angle(fourier)

    yy, xx = np.indices(a.shape)
    cy, cx = (ny - 1) / 2.0, (nx - 1) / 2.0
    kx = (xx - cx) / nx
    ky = (yy - cy) / ny
    kr = np.hypot(kx, ky)
    theta = np.arctan2(ky, kx)

    r_edges = np.linspace(0.02, 0.20, 7)
    a_edges = np.linspace(-math.pi, math.pi, 7)
    phases: list[float] = []
    amplitudes: list[float] = []
    counts: list[int] = []

    total_amp = float(np.sum(amp[(kr >= 0.02) & (kr < 0.20)]))
    for ri in range(6):
        for ai in range(6):
            mask = (
                (kr >= r_edges[ri])
                & (kr < r_edges[ri + 1])
                & (theta >= a_edges[ai])
                & (theta < a_edges[ai + 1])
            )
            w = amp[mask]
            p = phase[mask]
            counts.append(int(np.count_nonzero(mask)))
            if w.size == 0 or float(np.sum(w)) <= 1e-18:
                phases.append(0.0)
                amplitudes.append(0.0)
                continue
            z = np.sum(w * np.exp(1j * p))
            phases.append(float(np.angle(z) % (2.0 * math.pi)))
            amplitudes.append(float(abs(z) / np.sum(w)))

    if len(phases) != 36:
        raise RuntimeError("phase36 sectorizer did not produce exactly 36 coordinates")

    return {
        "phase36": phases,
        "sector_coherence36": amplitudes,
        "sector_counts36": counts,
        "global_sector_coherence": float(np.mean(amplitudes)),
        "fft_band_total_amplitude": total_amp,
        "layout": "RADIAL6_X_ANGULAR6_ROW_MAJOR",
        "radial_range_cycles_per_pixel": [0.02, 0.20],
        "authority": "MEASUREMENT_DERIVED_CANDIDATE",
        "semantic_axis_assignment": False,
    }


def circular_delta(left: float, right: float) -> float:
    return math.atan2(math.sin(left - right), math.cos(left - right))


def phase_rms_distance(left: list[float], right: list[float]) -> float:
    if len(left) != 36 or len(right) != 36:
        raise ValueError("GREMLIN phase distance requires exact 36D cells")
    return math.sqrt(
        sum(circular_delta(x, y) ** 2 for x, y in zip(left, right)) / 36.0
    )


def moire_overlap(left: list[float], right: list[float]) -> float:
    if len(left) != 36 or len(right) != 36:
        raise ValueError("GREMLIN moire overlap requires exact 36D cells")
    z = sum(
        np.exp(1j * circular_delta(x, y))
        for x, y in zip(left, right)
    ) / 36.0
    return float(abs(z))


def project_hologram(cells: list[list[float]]) -> dict[str, Any]:
    if not cells:
        raise ValueError("project_hologram requires at least one 36D cell")
    if any(len(cell) != 36 for cell in cells):
        raise ValueError("all hologram cells must be exact 36D")
    phases: list[float] = []
    amplitudes: list[float] = []
    for j in range(36):
        z = sum(np.exp(1j * cell[j]) for cell in cells)
        phases.append(float(np.angle(z) % (2.0 * math.pi)) if abs(z) > 1e-15 else 0.0)
        amplitudes.append(float(abs(z) / len(cells)))
    return {
        "phase36": phases,
        "amplitude36": amplitudes,
        "global_coherence": float(np.mean(amplitudes)),
        "source_count": len(cells),
        "authority": "CANDIDATE_ONLY",
    }


def center_stats(image: np.ndarray) -> dict[str, Any]:
    a = finite_image(image)
    ny, nx = a.shape
    yy, xx = np.indices(a.shape)
    cy, cx = (ny - 1) / 2.0, (nx - 1) / 2.0
    r_deg = np.hypot(xx - cx, yy - cy) * (FOV_DEG / WIDTH)

    aperture = r_deg <= 0.25
    annulus = (r_deg >= 0.75) & (r_deg <= 1.50)
    bg = a[annulus]
    bg_med = float(np.median(bg))
    mad = float(np.median(np.abs(bg - bg_med)))
    robust_sigma = 1.4826 * mad if mad > 0 else float(np.std(bg))
    center_mean = float(np.mean(a[aperture]))
    z = float((center_mean - bg_med) / robust_sigma) if robust_sigma > 0 else None

    return {
        "aperture_radius_deg": 0.25,
        "background_annulus_deg": [0.75, 1.50],
        "center_mean": center_mean,
        "background_median": bg_med,
        "background_robust_sigma": robust_sigma,
        "center_robust_z": z,
    }


def common_beam(image: np.ndarray, native_fwhm_arcmin: float) -> np.ndarray:
    if native_fwhm_arcmin >= COMMON_FWHM_ARCMIN:
        return finite_image(image)
    extra_fwhm = math.sqrt(
        COMMON_FWHM_ARCMIN**2 - native_fwhm_arcmin**2
    )
    sigma_pix = extra_fwhm / (2.354820045 * PIX_ARCMIN)
    a = finite_image(image)
    ny, nx = a.shape
    fy = np.fft.fftfreq(ny)
    fx = np.fft.fftfreq(nx)
    ky, kx = np.meshgrid(fy, fx, indexing="ij")
    transfer = np.exp(-2.0 * (math.pi**2) * (sigma_pix**2) * (kx*kx + ky*ky))
    return np.fft.ifft2(np.fft.fft2(a) * transfer).real


def corr(a: np.ndarray, b: np.ndarray) -> float | None:
    x = (finite_image(a) - np.mean(a)).ravel()
    y = (finite_image(b) - np.mean(b)).ravel()
    denom = float(np.linalg.norm(x) * np.linalg.norm(y))
    return float(np.dot(x, y) / denom) if denom > 0 else None


def fetch_cutout(source: dict[str, Any]) -> tuple[bytes, str]:
    params = {
        "hips": source["hips"],
        "width": str(WIDTH),
        "height": str(HEIGHT),
        "projection": "TAN",
        "fov": str(FOV_DEG),
        "ra": str(TARGET["ra_deg"]),
        "dec": str(TARGET["dec_deg"]),
        "coordsys": "icrs",
        "rotation_angle": "0",
        "format": "fits",
    }
    url = HIPS2FITS + "?" + urlencode(params)
    response = requests.get(
        url,
        timeout=120,
        headers={"User-Agent": "GREMLIN-FERRET-PBH-CMB/0.3 read-only"},
    )
    if response.status_code != 200:
        raise RuntimeError(f"HTTP {response.status_code}")
    if not response.content.startswith(b"SIMPLE"):
        head = response.content[:200].decode("latin1", "replace")
        raise RuntimeError(f"non-FITS response: {head}")
    return response.content, url


def main() -> int:
    receipt: dict[str, Any] = {
        "schema": "PBH_CMB_CENTERED_GREMLIN_T36_RECEIPT_V0_3",
        "target": TARGET,
        "cutout": {
            "width": WIDTH,
            "height": HEIGHT,
            "fov_deg": FOV_DEG,
            "pixel_scale_arcmin": PIX_ARCMIN,
            "projection": "TAN",
            "coordsys": "ICRS",
        },
        "common_beam_fwhm_arcmin": COMMON_FWHM_ARCMIN,
        "maps": {},
        "raw_pairwise_correlations": {},
        "common_beam_pairwise_correlations": {},
        "gremlin_source": GREMLIN_SOURCE,
        "gremlin_phase36_pairwise_common_beam": {},
        "gremlin_phase36_pairwise_raw": {},
        "gremlin_hologram_all_common_beam": None,
        "gremlin_hologram_frequency_common_beam": None,
        "status": "PARTIAL",
        "interpretation": "DESCRIPTIVE_ONLY_NOT_PBH_DETECTION",
        "authority": "CANDIDATE_ONLY",
        "canon_allowed": False,
    }

    raw_images: dict[str, np.ndarray] = {}
    smoothed_images: dict[str, np.ndarray] = {}
    phase36_raw: dict[str, list[float]] = {}
    phase36_common: dict[str, list[float]] = {}

    for name, source in SOURCES.items():
        entry: dict[str, Any] = {
            "hips_id": source["hips"],
            "kind": source["kind"],
            "native_fwhm_arcmin_approx": source["native_fwhm_arcmin"],
        }
        try:
            payload, url = fetch_cutout(source)
            entry["source_url"] = url
            entry["bytes"] = len(payload)
            entry["sha256"] = hashlib.sha256(payload).hexdigest()

            with fits.open(io.BytesIO(payload), memmap=False) as hdul:
                hdu = hdul[0]
                data = np.asarray(hdu.data, dtype=float)
                while data.ndim > 2:
                    data = data[0]
                if data.shape != (HEIGHT, WIDTH):
                    raise RuntimeError(f"unexpected cutout shape {data.shape}")
                data = finite_image(data)
                header = hdu.header

            smoothed = common_beam(data, source["native_fwhm_arcmin"])
            q = WIDTH // 4
            central = data[q:3*q, q:3*q]
            central_common = smoothed[q:3*q, q:3*q]

            entry.update(
                {
                    "header": {
                        key: (
                            list(header.get(key))
                            if key in ("HISTORY", "COMMENT")
                            else header.get(key)
                        )
                        for key in (
                            "BUNIT",
                            "CTYPE1",
                            "CTYPE2",
                            "CRVAL1",
                            "CRVAL2",
                            "CDELT1",
                            "CDELT2",
                            "NAXIS1",
                            "NAXIS2",
                        )
                        if key in header
                    },
                    "history": [str(v) for v in header.get("HISTORY", [])],
                    "center_stats_raw": center_stats(data),
                    "center_stats_common_beam": center_stats(smoothed),
                    "phase_markers_full_raw": phase_markers(data),
                    "phase_markers_central4deg_raw": phase_markers(central),
                    "phase_markers_central4deg_common_beam": phase_markers(
                        central_common
                    ),
                    "gremlin_phase36_central4deg_raw": phase36_sector_cell(
                        central
                    ),
                    "gremlin_phase36_central4deg_common_beam": phase36_sector_cell(
                        central_common
                    ),
                }
            )
            raw_images[name] = data
            smoothed_images[name] = smoothed
            phase36_raw[name] = entry["gremlin_phase36_central4deg_raw"]["phase36"]
            phase36_common[name] = entry[
                "gremlin_phase36_central4deg_common_beam"
            ]["phase36"]
            receipt["maps"][name] = entry
        except Exception as exc:
            entry["error"] = f"{type(exc).__name__}: {exc}"
            receipt["maps"][name] = entry

    names = sorted(raw_images)
    for i, left in enumerate(names):
        for right in names[i + 1 :]:
            key = f"{left}__{right}"
            receipt["raw_pairwise_correlations"][key] = corr(
                raw_images[left], raw_images[right]
            )
            receipt["common_beam_pairwise_correlations"][key] = corr(
                smoothed_images[left], smoothed_images[right]
            )
            receipt["gremlin_phase36_pairwise_raw"][key] = {
                "moire_overlap": moire_overlap(
                    phase36_raw[left], phase36_raw[right]
                ),
                "phase_rms_distance_rad": phase_rms_distance(
                    phase36_raw[left], phase36_raw[right]
                ),
            }
            receipt["gremlin_phase36_pairwise_common_beam"][key] = {
                "moire_overlap": moire_overlap(
                    phase36_common[left], phase36_common[right]
                ),
                "phase_rms_distance_rad": phase_rms_distance(
                    phase36_common[left], phase36_common[right]
                ),
            }

    if phase36_common:
        receipt["gremlin_hologram_all_common_beam"] = project_hologram(
            [phase36_common[name] for name in sorted(phase36_common)]
        )
        frequency_names = [
            name
            for name in sorted(phase36_common)
            if receipt["maps"][name].get("kind") == "CMB_FREQUENCY"
        ]
        if frequency_names:
            receipt["gremlin_hologram_frequency_common_beam"] = project_hologram(
                [phase36_common[name] for name in frequency_names]
            )

    ok = [k for k, v in receipt["maps"].items() if "error" not in v]
    receipt["successful_maps"] = ok
    receipt["status"] = "ANALYZED" if len(ok) == len(SOURCES) else "PARTIAL"

    payload = json.dumps(receipt, sort_keys=True, indent=2)
    with open("/tmp/receipt.json", "w", encoding="utf-8") as handle:
        handle.write(payload + "\n")

    compact = {
        "schema": receipt["schema"],
        "status": receipt["status"],
        "successful_maps": receipt["successful_maps"],
        "gremlin_hologram_all_common_beam": receipt[
            "gremlin_hologram_all_common_beam"
        ],
        "gremlin_hologram_frequency_common_beam": receipt[
            "gremlin_hologram_frequency_common_beam"
        ],
        "gremlin_phase36_pairwise_common_beam": receipt[
            "gremlin_phase36_pairwise_common_beam"
        ],
        "map_summary": {
            name: {
                "global_sector_coherence_common": row[
                    "gremlin_phase36_central4deg_common_beam"
                ]["global_sector_coherence"],
                "global_sector_coherence_raw": row[
                    "gremlin_phase36_central4deg_raw"
                ]["global_sector_coherence"],
                "center_robust_z_common": row[
                    "center_stats_common_beam"
                ]["center_robust_z"],
                "phase_coherence_common": row[
                    "phase_markers_central4deg_common_beam"
                ]["phase_coherence"],
                "phase_entropy_common": row[
                    "phase_markers_central4deg_common_beam"
                ]["phase_entropy"],
                "sha256": row["sha256"],
            }
            for name, row in receipt["maps"].items()
            if "error" not in row
        },
    }
    print(
        "PBH_CMB_T36_SUMMARY="
        + json.dumps(compact, sort_keys=True, separators=(",", ":")),
        flush=True,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
