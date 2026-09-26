#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import io
import json
import math
from concurrent.futures import ThreadPoolExecutor, as_completed
from itertools import combinations
from typing import Any
from urllib.parse import urlencode

import numpy as np
import requests
from astropy.coordinates import SkyCoord
from astropy.io import fits
import astropy.units as u

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
HIPS2FITS = "https://alasky.cds.unistra.fr/hips-image-services/hips2fits"

SOURCES = {
    "planck_100": {"hips": "CDS/P/PLANCK/R3/HFI100", "fwhm": 9.66},
    "planck_143": {"hips": "CDS/P/PLANCK/R3/HFI143", "fwhm": 7.27},
    "planck_217": {"hips": "CDS/P/PLANCK/R3/HFI217", "fwhm": 5.01},
    "wmap_V": {"hips": "CDS/P/WMAP/V/9yr", "fwhm": 21.0},
    "wmap_W": {"hips": "CDS/P/WMAP/W/9yr", "fwhm": 13.2},
}

CONTROL_LONGITUDES = [
    0.0, 15.0, 30.0, 45.0, 60.0, 75.0, 90.0, 105.0, 120.0,
    135.0, 150.0, 165.0, 180.0, 195.0, 270.0, 285.0, 300.0,
    315.0, 330.0, 345.0,
]

PRIMARY_STATISTIC = "median_pairwise_moire_frequency_maps_common_beam"
INTERPRETATION = (
    "EXPLORATORY_POSTHOC_MATCHED_LATITUDE_NULL_ONLY_NOT_FORMAL_PVALUE"
)


def finite_image(image: np.ndarray) -> np.ndarray:
    image = np.asarray(image, dtype=float)
    finite = np.isfinite(image)
    if not finite.any():
        raise RuntimeError("image contains no finite pixels")
    med = float(np.nanmedian(image[finite]))
    return np.where(finite, image, med)


def common_beam(image: np.ndarray, native_fwhm_arcmin: float) -> np.ndarray:
    a = finite_image(image)
    if native_fwhm_arcmin >= COMMON_FWHM_ARCMIN:
        return a
    extra_fwhm = math.sqrt(COMMON_FWHM_ARCMIN**2 - native_fwhm_arcmin**2)
    sigma_pix = extra_fwhm / (2.354820045 * PIX_ARCMIN)
    ny, nx = a.shape
    fy = np.fft.fftfreq(ny)
    fx = np.fft.fftfreq(nx)
    ky, kx = np.meshgrid(fy, fx, indexing="ij")
    transfer = np.exp(
        -2.0 * (math.pi**2) * (sigma_pix**2) * (kx*kx + ky*ky)
    )
    return np.fft.ifft2(np.fft.fft2(a) * transfer).real


def phase36_sector_cell(image: np.ndarray) -> list[float]:
    a = finite_image(image)
    centered = a - float(np.median(a))
    ny, nx = a.shape
    window = np.outer(np.hanning(ny), np.hanning(nx))
    F = np.fft.fftshift(np.fft.fft2(centered * window))
    amp = np.abs(F)
    phase = np.angle(F)

    yy, xx = np.indices(a.shape)
    cy, cx = (ny - 1) / 2.0, (nx - 1) / 2.0
    kx = (xx - cx) / nx
    ky = (yy - cy) / ny
    kr = np.hypot(kx, ky)
    theta = np.arctan2(ky, kx)

    r_edges = np.linspace(0.02, 0.20, 7)
    a_edges = np.linspace(-math.pi, math.pi, 7)
    phases: list[float] = []
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
            if w.size == 0 or float(np.sum(w)) <= 1e-18:
                phases.append(0.0)
                continue
            z = np.sum(w * np.exp(1j * p))
            phases.append(float(np.angle(z) % (2.0 * math.pi)))
    if len(phases) != 36:
        raise RuntimeError("T36 sectorizer did not produce exactly 36 values")
    return phases


def circular_delta(left: float, right: float) -> float:
    return math.atan2(math.sin(left - right), math.cos(left - right))


def moire_overlap(left: list[float], right: list[float]) -> float:
    if len(left) != 36 or len(right) != 36:
        raise ValueError("moire_overlap requires exact T36 cells")
    z = sum(
        np.exp(1j * circular_delta(x, y))
        for x, y in zip(left, right)
    ) / 36.0
    return float(abs(z))


def fetch_map(name: str, source: dict[str, Any], ra: float, dec: float) -> tuple[str, np.ndarray, dict[str, Any]]:
    params = {
        "hips": source["hips"],
        "width": str(WIDTH),
        "height": str(HEIGHT),
        "projection": "TAN",
        "fov": str(FOV_DEG),
        "ra": f"{ra:.12f}",
        "dec": f"{dec:.12f}",
        "coordsys": "icrs",
        "rotation_angle": "0",
        "format": "fits",
    }
    url = HIPS2FITS + "?" + urlencode(params)
    response = requests.get(
        url,
        timeout=30,
        headers={"User-Agent": "GREMLIN-FERRET-PBH-CMB-NULL/0.1 read-only"},
    )
    if response.status_code != 200 or not response.content.startswith(b"SIMPLE"):
        raise RuntimeError(
            f"{name}: HTTP={response.status_code} bytes={len(response.content)}"
        )
    with fits.open(io.BytesIO(response.content), memmap=False) as hdul:
        data = np.asarray(hdul[0].data, dtype=float)
        while data.ndim > 2:
            data = data[0]
    if data.shape != (HEIGHT, WIDTH):
        raise RuntimeError(f"{name}: unexpected shape {data.shape}")
    meta = {
        "sha256": hashlib.sha256(response.content).hexdigest(),
        "bytes": len(response.content),
        "source_url": url,
    }
    return name, finite_image(data), meta


def analyze_coordinate(ra: float, dec: float) -> dict[str, Any]:
    raw: dict[str, np.ndarray] = {}
    meta: dict[str, Any] = {}
    with ThreadPoolExecutor(max_workers=5) as pool:
        futs = {
            pool.submit(fetch_map, name, source, ra, dec): name
            for name, source in SOURCES.items()
        }
        for fut in as_completed(futs):
            name, image, row = fut.result()
            raw[name] = image
            meta[name] = row

    cells: dict[str, list[float]] = {}
    for name, source in SOURCES.items():
        smooth = common_beam(raw[name], source["fwhm"])
        q = WIDTH // 4
        central = smooth[q:3*q, q:3*q]
        cells[name] = phase36_sector_cell(central)

    overlaps: dict[str, float] = {}
    values: list[float] = []
    for left, right in combinations(sorted(cells), 2):
        value = moire_overlap(cells[left], cells[right])
        overlaps[f"{left}__{right}"] = value
        values.append(value)

    arr = np.asarray(values, dtype=float)
    return {
        "statistic": {
            "median_pairwise_moire": float(np.median(arr)),
            "mean_pairwise_moire": float(np.mean(arr)),
            "min_pairwise_moire": float(np.min(arr)),
            "max_pairwise_moire": float(np.max(arr)),
            "std_pairwise_moire": float(np.std(arr)),
        },
        "pairwise_moire": overlaps,
        "hashes": {name: meta[name]["sha256"] for name in sorted(meta)},
    }


def galactic_to_icrs(l_deg: float, b_deg: float) -> tuple[float, float]:
    c = SkyCoord(l=l_deg*u.deg, b=b_deg*u.deg, frame="galactic").icrs
    return float(c.ra.deg), float(c.dec.deg)


def angular_lon_sep(a: float, b: float) -> float:
    d = abs((a - b + 180.0) % 360.0 - 180.0)
    return d


def empirical_upper_rank(target: float, controls: list[float]) -> float:
    # Add-one rank; descriptive only because the statistic was selected after
    # inspecting the target in an earlier exploratory pass.
    return float((1 + sum(v >= target for v in controls)) / (len(controls) + 1))


def main() -> int:
    print("T36_NULL_START target", TARGET, flush=True)
    target = analyze_coordinate(TARGET["ra_deg"], TARGET["dec_deg"])
    print(
        "T36_NULL_TARGET="
        + json.dumps(target["statistic"], sort_keys=True, separators=(",", ":")),
        flush=True,
    )

    controls = []
    for i, l_deg in enumerate(CONTROL_LONGITUDES, 1):
        if angular_lon_sep(l_deg, TARGET["gal_l_deg"]) < 30.0:
            continue
        ra, dec = galactic_to_icrs(l_deg, TARGET["gal_b_deg"])
        row = analyze_coordinate(ra, dec)
        controls.append(
            {
                "l_deg": l_deg,
                "b_deg": TARGET["gal_b_deg"],
                "ra_deg": ra,
                "dec_deg": dec,
                **row["statistic"],
            }
        )
        print(
            f"T36_NULL_CONTROL {len(controls)} l={l_deg} median={row['statistic']['median_pairwise_moire']:.9f}",
            flush=True,
        )

    medians = [r["median_pairwise_moire"] for r in controls]
    means = [r["mean_pairwise_moire"] for r in controls]
    minima = [r["min_pairwise_moire"] for r in controls]
    receipt = {
        "schema": "PBH_CMB_T36_MATCHED_LATITUDE_NULL_V0_1",
        "status": "ANALYZED",
        "target": {
            **TARGET,
            **target["statistic"],
            "pairwise_moire": target["pairwise_moire"],
            "hashes": target["hashes"],
        },
        "control_definition": {
            "same_galactic_latitude_deg": TARGET["gal_b_deg"],
            "control_longitudes_deg": CONTROL_LONGITUDES,
            "target_exclusion_deg": 30.0,
            "common_beam_fwhm_arcmin": COMMON_FWHM_ARCMIN,
            "maps": sorted(SOURCES),
            "t36_layout": "RADIAL6_X_ANGULAR6_ROW_MAJOR",
        },
        "n_controls": len(controls),
        "controls": controls,
        "empirical_descriptive_ranks": {
            "median_pairwise_moire_upper": empirical_upper_rank(
                target["statistic"]["median_pairwise_moire"], medians
            ),
            "mean_pairwise_moire_upper": empirical_upper_rank(
                target["statistic"]["mean_pairwise_moire"], means
            ),
            "min_pairwise_moire_upper": empirical_upper_rank(
                target["statistic"]["min_pairwise_moire"], minima
            ),
        },
        "primary_statistic": PRIMARY_STATISTIC,
        "interpretation": INTERPRETATION,
        "authority": "CANDIDATE_ONLY",
        "canon_allowed": False,
    }

    compact = {
        "schema": receipt["schema"],
        "status": receipt["status"],
        "n_controls": receipt["n_controls"],
        "target_statistic": target["statistic"],
        "empirical_descriptive_ranks": receipt["empirical_descriptive_ranks"],
        "control_median_pairwise_moire": {
            "min": float(np.min(medians)),
            "median": float(np.median(medians)),
            "max": float(np.max(medians)),
            "mean": float(np.mean(medians)),
        },
        "interpretation": INTERPRETATION,
    }
    print(
        "PBH_T36_NULL_SUMMARY="
        + json.dumps(compact, sort_keys=True, separators=(",", ":")),
        flush=True,
    )
    print(
        "PBH_T36_NULL_RECEIPT="
        + json.dumps(receipt, sort_keys=True, separators=(",", ":")),
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
