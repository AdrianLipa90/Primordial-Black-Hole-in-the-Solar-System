#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import io
import json
import math
import sys
from itertools import combinations
from typing import Any
from urllib.parse import urlencode

import numpy as np
import requests
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.io import fits

WIDTH = 192
HEIGHT = 192
FOV_DEG = 8.0
PIX_ARCMIN = FOV_DEG * 60.0 / WIDTH
COMMON_FWHM_ARCMIN = 40.0
TARGET_GAL_L = 238.1709218507047
TARGET_GAL_B = -7.665782897020605

CONTROL_L = [0.0, 45.0, 90.0, 135.0, 180.0, 270.0, 315.0]

SOURCES = {
    "planck_smica": {"hips": "CDS/P/PLANCK/R3/CMB", "kind": "CMB_COMPONENT", "native_fwhm_arcmin": 5.0},
    "planck_100": {"hips": "CDS/P/PLANCK/R3/HFI100", "kind": "CMB_FREQUENCY", "native_fwhm_arcmin": 9.66},
    "planck_143": {"hips": "CDS/P/PLANCK/R3/HFI143", "kind": "CMB_FREQUENCY", "native_fwhm_arcmin": 7.27},
    "planck_217": {"hips": "CDS/P/PLANCK/R3/HFI217", "kind": "CMB_FREQUENCY", "native_fwhm_arcmin": 5.01},
    "wmap_V": {"hips": "CDS/P/WMAP/V/9yr", "kind": "CMB_FREQUENCY", "native_fwhm_arcmin": 21.0},
    "wmap_W": {"hips": "CDS/P/WMAP/W/9yr", "kind": "CMB_FREQUENCY", "native_fwhm_arcmin": 13.2},
}
HIPS2FITS = "https://alasky.cds.unistra.fr/hips-image-services/hips2fits"


def finite_image(image: np.ndarray) -> np.ndarray:
    a = np.asarray(image, dtype=float)
    m = np.isfinite(a)
    if not m.any():
        raise RuntimeError("no finite pixels")
    med = float(np.nanmedian(a[m]))
    return np.where(m, a, med)


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
    transfer = np.exp(-2.0 * math.pi**2 * sigma_pix**2 * (kx*kx + ky*ky))
    return np.fft.ifft2(np.fft.fft2(a) * transfer).real


def phase36_sector_cell(image: np.ndarray) -> list[float]:
    a = finite_image(image)
    ny, nx = a.shape
    centered = a - float(np.median(a))
    F = np.fft.fftshift(np.fft.fft2(centered * np.outer(np.hanning(ny), np.hanning(nx))))
    amp = np.abs(F)
    ph = np.angle(F)
    yy, xx = np.indices(a.shape)
    cy, cx = (ny - 1) / 2.0, (nx - 1) / 2.0
    kx = (xx - cx) / nx
    ky = (yy - cy) / ny
    kr = np.hypot(kx, ky)
    th = np.arctan2(ky, kx)
    r_edges = np.linspace(0.02, 0.20, 7)
    a_edges = np.linspace(-math.pi, math.pi, 7)
    out = []
    for ri in range(6):
        for ai in range(6):
            mask = (
                (kr >= r_edges[ri]) & (kr < r_edges[ri + 1])
                & (th >= a_edges[ai]) & (th < a_edges[ai + 1])
            )
            w = amp[mask]
            p = ph[mask]
            if w.size == 0 or float(np.sum(w)) <= 1e-18:
                out.append(0.0)
            else:
                z = np.sum(w * np.exp(1j * p))
                out.append(float(np.angle(z) % (2.0 * math.pi)))
    if len(out) != 36:
        raise RuntimeError("T36 sectorizer shape failure")
    return out


def circular_delta(a: float, b: float) -> float:
    return math.atan2(math.sin(a-b), math.cos(a-b))


def moire_overlap(left: list[float], right: list[float]) -> float:
    z = sum(np.exp(1j * circular_delta(a,b)) for a,b in zip(left,right)) / 36.0
    return float(abs(z))


def project_hologram(cells: list[list[float]]) -> float:
    amps = []
    for j in range(36):
        z = sum(np.exp(1j * cell[j]) for cell in cells)
        amps.append(abs(z) / len(cells))
    return float(np.mean(amps))


def corr(a: np.ndarray, b: np.ndarray) -> float | None:
    x=(finite_image(a)-np.mean(a)).ravel()
    y=(finite_image(b)-np.mean(b)).ravel()
    d=float(np.linalg.norm(x)*np.linalg.norm(y))
    return float(np.dot(x,y)/d) if d>0 else None


def center_z(image: np.ndarray) -> float | None:
    a=finite_image(image)
    ny,nx=a.shape
    yy,xx=np.indices(a.shape)
    cy,cx=(ny-1)/2.0,(nx-1)/2.0
    rdeg=np.hypot(xx-cx,yy-cy)*(FOV_DEG/WIDTH)
    ap=rdeg<=0.25
    ann=(rdeg>=0.75)&(rdeg<=1.50)
    bg=a[ann]
    med=float(np.median(bg))
    mad=float(np.median(np.abs(bg-med)))
    sig=1.4826*mad if mad>0 else float(np.std(bg))
    return float((np.mean(a[ap])-med)/sig) if sig>0 else None


def gal_to_icrs(l_deg: float, b_deg: float) -> tuple[float,float]:
    c=SkyCoord(l=l_deg*u.deg,b=b_deg*u.deg,frame="galactic").icrs
    return float(c.ra.deg),float(c.dec.deg)


def fetch(source: dict[str,Any], ra: float, dec: float) -> tuple[np.ndarray,str,str]:
    params={
        "hips":source["hips"],"width":str(WIDTH),"height":str(HEIGHT),
        "projection":"TAN","fov":str(FOV_DEG),"ra":str(ra),"dec":str(dec),
        "coordsys":"icrs","rotation_angle":"0","format":"fits"
    }
    url=HIPS2FITS+"?"+urlencode(params)
    r=requests.get(url,timeout=120,headers={"User-Agent":"GREMLIN-FERRET-PBH-CMB-null/0.2 read-only"})
    if r.status_code!=200 or not r.content.startswith(b"SIMPLE"):
        raise RuntimeError(f"source failure {source['hips']} HTTP {r.status_code}")
    with fits.open(io.BytesIO(r.content),memmap=False) as hdul:
        data=np.asarray(hdul[0].data,dtype=float)
        while data.ndim>2: data=data[0]
    if data.shape!=(HEIGHT,WIDTH):
        raise RuntimeError(f"unexpected shape {data.shape}")
    return finite_image(data),hashlib.sha256(r.content).hexdigest(),url


def analyze_point(l_deg: float, b_deg: float, label: str) -> dict[str,Any]:
    ra,dec=gal_to_icrs(l_deg,b_deg)
    images={}
    cells={}
    hashes={}
    urls={}
    for name,src in SOURCES.items():
        raw,sha,url=fetch(src,ra,dec)
        sm=common_beam(raw,src["native_fwhm_arcmin"])
        q=WIDTH//4
        central=sm[q:3*q,q:3*q]
        images[name]=sm
        cells[name]=phase36_sector_cell(central)
        hashes[name]=sha
        urls[name]=url

    pair={}
    vals=[]
    for left,right in combinations(sorted(cells),2):
        ov=moire_overlap(cells[left],cells[right])
        pair[f"{left}__{right}"]=ov
        vals.append(ov)

    freq=[cells[name] for name in sorted(cells) if SOURCES[name]["kind"]=="CMB_FREQUENCY"]
    all_cells=[cells[name] for name in sorted(cells)]

    return {
        "label":label,
        "galactic":{"l_deg":l_deg,"b_deg":b_deg},
        "icrs":{"ra_deg":ra,"dec_deg":dec},
        "mean_pairwise_moire":float(np.mean(vals)),
        "median_pairwise_moire":float(np.median(vals)),
        "hologram_all_global_coherence":project_hologram(all_cells),
        "hologram_frequency_global_coherence":project_hologram(freq),
        "smica_center_z":center_z(images["planck_smica"]),
        "smica_wmapV_corr":corr(images["planck_smica"],images["wmap_V"]),
        "smica_wmapW_corr":corr(images["planck_smica"],images["wmap_W"]),
        "wmapVW_corr":corr(images["wmap_V"],images["wmap_W"]),
        "pairwise_moire":pair,
        "hashes":hashes,
        "source_urls":urls,
    }


def empirical_upper(target: float, controls: list[float]) -> float:
    return (1.0 + sum(v >= target for v in controls)) / (len(controls)+1.0)


def empirical_abs_upper(target: float, controls: list[float]) -> float:
    return (1.0 + sum(abs(v) >= abs(target) for v in controls)) / (len(controls)+1.0)


def main() -> int:
    target=analyze_point(TARGET_GAL_L,TARGET_GAL_B,"PBH_TARGET")
    controls=[]
    for i,l in enumerate(CONTROL_L,1):
        print(f"CONTROL {i}/{len(CONTROL_L)} l={l}",flush=True)
        controls.append(analyze_point(l,TARGET_GAL_B,f"CONTROL_L{l:g}"))

    receipt={
        "schema":"PBH_CMB_GREMLIN_T36_MATCHED_LATITUDE_NULL_V0_2",
        "status":"ANALYZED",
        "target":target,
        "controls":controls,
        "control_definition":{
            "same_galactic_latitude_deg":TARGET_GAL_B,
            "longitude_grid_deg":CONTROL_L,
            "target_longitude_deg":TARGET_GAL_L,
            "field_deg":FOV_DEG,
            "common_beam_fwhm_arcmin":COMMON_FWHM_ARCMIN,
            "t36_layout":"RADIAL6_X_ANGULAR6",
        },
        "empirical_upper_tail_plus1":{
            "mean_pairwise_moire":empirical_upper(target["mean_pairwise_moire"],[x["mean_pairwise_moire"] for x in controls]),
            "hologram_all_global_coherence":empirical_upper(target["hologram_all_global_coherence"],[x["hologram_all_global_coherence"] for x in controls]),
            "hologram_frequency_global_coherence":empirical_upper(target["hologram_frequency_global_coherence"],[x["hologram_frequency_global_coherence"] for x in controls]),
            "abs_smica_center_z":empirical_abs_upper(target["smica_center_z"],[x["smica_center_z"] for x in controls]),
            "smica_wmapV_corr":empirical_upper(target["smica_wmapV_corr"],[x["smica_wmapV_corr"] for x in controls]),
            "smica_wmapW_corr":empirical_upper(target["smica_wmapW_corr"],[x["smica_wmapW_corr"] for x in controls]),
        },
        "interpretation":"PILOT_MATCHED_LATITUDE_NULL_ONLY_NOT_FORMAL_SURVEY_PVALUE_NOT_PBH_DETECTION",
        "gremlin_algebra_source":{
            "repository":"AdrianLipa90/GREMLIN",
            "commit":"9fdc42e507d5fe54c13dda975c297173ba0e4eb6",
            "artifact":"gremlin_mcp/holographic_orbital_memory.py",
            "operators":["moire_overlap","project_hologram"],
            "phase36_mapping":"measurement-derived fixed 6x6 Fourier sector cell"
        },
        "authority":"CANDIDATE_ONLY",
        "canon_allowed":False,
    }
    print("PBH_T36_NULL_RECEIPT="+json.dumps(receipt,sort_keys=True,separators=(",",":")),flush=True)
    return 0


if __name__=="__main__":
    sys.exit(main())
