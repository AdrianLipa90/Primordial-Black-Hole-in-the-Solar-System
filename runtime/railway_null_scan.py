#!/usr/bin/env python3
from __future__ import annotations

import io
import json
import math
import sys
from urllib.parse import urlencode

import numpy as np
import requests
from astropy.coordinates import SkyCoord
from astropy.io import fits
import astropy.units as u

TARGET_L = 238.1709218507047
TARGET_B = -7.665782897020605
WIDTH = 256
HEIGHT = 256
FOV_DEG = 8.0
PIX_ARCMIN = FOV_DEG * 60.0 / WIDTH
COMMON_FWHM_ARCMIN = 40.0
HIPS2FITS = "https://alasky.cds.unistra.fr/hips-image-services/hips2fits"

SOURCES = {
    "planck_smica": {
        "hips": "CDS/P/PLANCK/R3/CMB",
        "native_fwhm_arcmin": 5.0,
    },
    "planck_143": {
        "hips": "CDS/P/PLANCK/R3/HFI143",
        "native_fwhm_arcmin": 7.27,
    },
    "wmap_W": {
        "hips": "CDS/P/WMAP/W/9yr",
        "native_fwhm_arcmin": 13.2,
    },
}

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


def finite_image(image: np.ndarray) -> np.ndarray:
    a = np.asarray(image, dtype=float)
    finite = np.isfinite(a)
    if not finite.any():
        raise RuntimeError("no finite pixels")
    med = float(np.nanmedian(a[finite]))
    return np.where(finite, a, med)


def common_beam(image: np.ndarray, native_fwhm_arcmin: float) -> np.ndarray:
    a = finite_image(image)
    if native_fwhm_arcmin >= COMMON_FWHM_ARCMIN:
        return a
    extra_fwhm = math.sqrt(
        COMMON_FWHM_ARCMIN**2 - native_fwhm_arcmin**2
    )
    sigma_pix = extra_fwhm / (2.354820045 * PIX_ARCMIN)
    ny, nx = a.shape
    fy = np.fft.fftfreq(ny)
    fx = np.fft.fftfreq(nx)
    ky, kx = np.meshgrid(fy, fx, indexing="ij")
    transfer = np.exp(
        -2.0 * (math.pi**2) * (sigma_pix**2) * (kx * kx + ky * ky)
    )
    return np.fft.ifft2(np.fft.fft2(a) * transfer).real


def phase36_sector_cell(image: np.ndarray) -> dict:
    """Measurement-derived exact T^36: 6 radial x 6 angular Fourier sectors."""
    a = finite_image(image)
    ny, nx = a.shape
    centered = a - float(np.median(a))
    fourier = np.fft.fftshift(
        np.fft.fft2(centered * np.outer(np.hanning(ny), np.hanning(nx)))
    )
    amp = np.abs(fourier)
    phase = np.angle(fourier)

    yy, xx = np.indices(a.shape)
    cy, cx = (ny - 1) / 2.0, (nx - 1) / 2.0
    kx = (xx - cx) / nx
    ky = (yy - cy) / ny
    kr = np.hypot(kx, ky)
    theta = np.arctan2(ky, kx)

    radial_edges = np.linspace(0.02, 0.20, 7)
    angle_edges = np.linspace(-math.pi, math.pi, 7)
    phases: list[float] = []
    coherences: list[float] = []

    for ri in range(6):
        for ai in range(6):
            mask = (
                (kr >= radial_edges[ri])
                & (kr < radial_edges[ri + 1])
                & (theta >= angle_edges[ai])
                & (theta < angle_edges[ai + 1])
            )
            w = amp[mask]
            p = phase[mask]
            if w.size == 0 or float(np.sum(w)) <= 1e-18:
                phases.append(0.0)
                coherences.append(0.0)
                continue
            z = np.sum(w * np.exp(1j * p))
            phases.append(float(np.angle(z) % (2.0 * math.pi)))
            coherences.append(float(abs(z) / np.sum(w)))

    if len(phases) != 36:
        raise RuntimeError("phase36 sectorizer did not produce exactly 36 coordinates")

    return {
        "phase36": phases,
        "global_sector_coherence": float(np.mean(coherences)),
        "semantic_axis_assignment": False,
        "layout": "RADIAL6_X_ANGULAR6_ROW_MAJOR",
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


def project_hologram(cells: list[list[float]]) -> dict:
    if not cells or any(len(cell) != 36 for cell in cells):
        raise ValueError("project_hologram requires exact 36D cells")
    phases: list[float] = []
    amplitudes: list[float] = []
    for j in range(36):
        z = sum(np.exp(1j * cell[j]) for cell in cells)
        phases.append(
            float(np.angle(z) % (2.0 * math.pi)) if abs(z) > 1e-15 else 0.0
        )
        amplitudes.append(float(abs(z) / len(cells)))
    return {
        "phase36": phases,
        "amplitude36": amplitudes,
        "global_coherence": float(np.mean(amplitudes)),
        "source_count": len(cells),
    }


def corr(a: np.ndarray, b: np.ndarray) -> float | None:
    x = (finite_image(a) - float(np.mean(a))).ravel()
    y = (finite_image(b) - float(np.mean(b))).ravel()
    denom = float(np.linalg.norm(x) * np.linalg.norm(y))
    return float(np.dot(x, y) / denom) if denom > 0 else None


def center_stats(image: np.ndarray) -> dict:
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
    center_z = (
        float((center_mean - bg_med) / robust_sigma)
        if robust_sigma > 0
        else None
    )
    return {
        "center_mean": center_mean,
        "background_median": bg_med,
        "background_robust_sigma": robust_sigma,
        "center_robust_z": center_z,
    }


def gal_to_icrs(l_deg: float, b_deg: float) -> tuple[float, float]:
    c = SkyCoord(l=l_deg * u.deg, b=b_deg * u.deg, frame="galactic").icrs
    return float(c.ra.deg), float(c.dec.deg)


def fetch_cutout(hips: str, ra_deg: float, dec_deg: float) -> np.ndarray:
    params = {
        "hips": hips,
        "width": WIDTH,
        "height": HEIGHT,
        "projection": "TAN",
        "fov": FOV_DEG,
        "ra": ra_deg,
        "dec": dec_deg,
        "coordsys": "icrs",
        "rotation_angle": 0,
        "format": "fits",
    }
    url = HIPS2FITS + "?" + urlencode(params)
    response = requests.get(
        url,
        timeout=90,
        headers={"User-Agent": "GREMLIN-FERRET-PBH-NULL/0.2 read-only"},
    )
    if response.status_code != 200 or not response.content.startswith(b"SIMPLE"):
        raise RuntimeError(
            f"fetch failed status={response.status_code} "
            f"head={response.content[:120]!r}"
        )
    with fits.open(io.BytesIO(response.content), memmap=False) as hdul:
        data = np.asarray(hdul[0].data, dtype=float)
        while data.ndim > 2:
            data = data[0]
    if data.shape != (HEIGHT, WIDTH):
        raise RuntimeError(f"unexpected shape {data.shape}")
    return finite_image(data)


def evaluate_position(label: str, l_deg: float, b_deg: float) -> dict:
    ra_deg, dec_deg = gal_to_icrs(l_deg, b_deg)
    maps: dict[str, dict] = {}
    cells: dict[str, list[float]] = {}
    images: dict[str, np.ndarray] = {}

    for name, source in SOURCES.items():
        raw = fetch_cutout(source["hips"], ra_deg, dec_deg)
        smoothed = common_beam(raw, source["native_fwhm_arcmin"])
        q = WIDTH // 4
        central = smoothed[q:3*q, q:3*q]
        t36 = phase36_sector_cell(central)
        maps[name] = {
            "center_stats": center_stats(smoothed),
            "phase36_global_sector_coherence": t36["global_sector_coherence"],
        }
        cells[name] = t36["phase36"]
        images[name] = central

    names = sorted(cells)
    pairwise_moire: dict[str, float] = {}
    pairwise_phase_rms: dict[str, float] = {}
    pairwise_pixel_corr: dict[str, float | None] = {}

    for i, left in enumerate(names):
        for right in names[i + 1 :]:
            key = f"{left}__{right}"
            pairwise_moire[key] = moire_overlap(cells[left], cells[right])
            pairwise_phase_rms[key] = phase_rms_distance(
                cells[left], cells[right]
            )
            pairwise_pixel_corr[key] = corr(images[left], images[right])

    hologram = project_hologram([cells[name] for name in names])

    return {
        "label": label,
        "galactic": {"l_deg": l_deg % 360.0, "b_deg": b_deg},
        "icrs": {"ra_deg": ra_deg, "dec_deg": dec_deg},
        "maps": maps,
        "pairwise_moire": pairwise_moire,
        "pairwise_phase_rms_rad": pairwise_phase_rms,
        "pairwise_pixel_corr": pairwise_pixel_corr,
        "mean_pairwise_moire": float(np.mean(list(pairwise_moire.values()))),
        "mean_pairwise_pixel_corr": float(
            np.mean([v for v in pairwise_pixel_corr.values() if v is not None])
        ),
        "hologram_global_coherence": hologram["global_coherence"],
        "smica_abs_center_z": abs(
            maps["planck_smica"]["center_stats"]["center_robust_z"]
        ),
    }


def empirical_upper_tail(target: float, controls: list[float]) -> dict:
    ge = sum(1 for value in controls if value >= target)
    return {
        "target": target,
        "controls_ge_target": ge,
        "control_count": len(controls),
        "rank_from_top": ge + 1,
        "empirical_upper_tail_fraction": (ge + 1) / (len(controls) + 1),
    }


def main() -> int:
    positions = [("TARGET", TARGET_L, TARGET_B)]
    for k in range(1, 12):
        positions.append(
            (f"CTRL_{k:02d}", (TARGET_L + 30.0 * k) % 360.0, TARGET_B)
        )

    rows: list[dict] = []
    for index, (label, l_deg, b_deg) in enumerate(positions):
        print(
            f"NULLSCAN {index + 1}/{len(positions)} "
            f"{label} l={l_deg:.6f} b={b_deg:.6f}",
            flush=True,
        )
        rows.append(evaluate_position(label, l_deg, b_deg))

    target = rows[0]
    controls = rows[1:]

    metrics = {
        "mean_pairwise_moire": empirical_upper_tail(
            target["mean_pairwise_moire"],
            [row["mean_pairwise_moire"] for row in controls],
        ),
        "hologram_global_coherence": empirical_upper_tail(
            target["hologram_global_coherence"],
            [row["hologram_global_coherence"] for row in controls],
        ),
        "mean_pairwise_pixel_corr": empirical_upper_tail(
            target["mean_pairwise_pixel_corr"],
            [row["mean_pairwise_pixel_corr"] for row in controls],
        ),
        "smica_abs_center_z": empirical_upper_tail(
            target["smica_abs_center_z"],
            [row["smica_abs_center_z"] for row in controls],
        ),
    }

    receipt = {
        "schema": "PBH_CMB_MATCHED_GALACTIC_NULL_SCAN_V0_1",
        "status": "ANALYZED",
        "target": {"l_deg": TARGET_L, "b_deg": TARGET_B},
        "control_design": (
            "11 deterministic controls at identical Galactic latitude; "
            "longitude = target + 30*k degrees, k=1..11"
        ),
        "maps": list(SOURCES),
        "cutout": {
            "fov_deg": FOV_DEG,
            "width": WIDTH,
            "height": HEIGHT,
            "common_beam_fwhm_arcmin": COMMON_FWHM_ARCMIN,
        },
        "gremlin_source": GREMLIN_SOURCE,
        "target_row": target,
        "control_rows": controls,
        "target_empirical_ranks": metrics,
        "interpretation": (
            "MATCHED_CONTROL_CALIBRATION_ONLY_NOT_PBH_DETECTION"
        ),
        "authority": "CANDIDATE_ONLY",
        "canon_allowed": False,
    }

    print(
        "PBH_NULL_RECEIPT="
        + json.dumps(receipt, sort_keys=True, separators=(",", ":")),
        flush=True,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
