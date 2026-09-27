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
from astropy.coordinates import SkyCoord
from astropy.io import fits
import astropy.units as u

TARGET = {
    "ra_deg": 107.576646,
    "dec_deg": -25.907579,
    "gal_l_deg": 238.1709218507047,
    "gal_b_deg": -7.665782897020605,
}

WIDTH = 128
HEIGHT = 128
FOV_DEG = 6.0
PIX_ARCMIN = FOV_DEG * 60.0 / WIDTH
COMMON_FWHM_ARCMIN = 40.0
EXCLUSION_DEG = 30.0

SOURCES = {
    "planck_smica": {"hips":"CDS/P/PLANCK/R3/CMB","kind":"CMB_COMPONENT","native_fwhm_arcmin":5.0},
    "planck_100": {"hips":"CDS/P/PLANCK/R3/HFI100","kind":"CMB_FREQUENCY","native_fwhm_arcmin":9.66},
    "planck_143": {"hips":"CDS/P/PLANCK/R3/HFI143","kind":"CMB_FREQUENCY","native_fwhm_arcmin":7.27},
    "planck_217": {"hips":"CDS/P/PLANCK/R3/HFI217","kind":"CMB_FREQUENCY","native_fwhm_arcmin":5.01},
    "wmap_V": {"hips":"CDS/P/WMAP/V/9yr","kind":"CMB_FREQUENCY","native_fwhm_arcmin":21.0},
    "wmap_W": {"hips":"CDS/P/WMAP/W/9yr","kind":"CMB_FREQUENCY","native_fwhm_arcmin":13.2},
}
HIPS2FITS = "https://alasky.cds.unistra.fr/hips-image-services/hips2fits"


def finite_image(image: np.ndarray) -> np.ndarray:
    a=np.asarray(image,dtype=float)
    m=np.isfinite(a)
    if not m.any():
        raise RuntimeError("no finite pixels")
    med=float(np.nanmedian(a[m]))
    return np.where(m,a,med)


def common_beam(image: np.ndarray, native_fwhm_arcmin: float) -> np.ndarray:
    a=finite_image(image)
    if native_fwhm_arcmin >= COMMON_FWHM_ARCMIN:
        return a
    extra=math.sqrt(COMMON_FWHM_ARCMIN**2-native_fwhm_arcmin**2)
    sigma_pix=extra/(2.354820045*PIX_ARCMIN)
    ny,nx=a.shape
    fy=np.fft.fftfreq(ny); fx=np.fft.fftfreq(nx)
    ky,kx=np.meshgrid(fy,fx,indexing="ij")
    transfer=np.exp(-2.0*(math.pi**2)*(sigma_pix**2)*(kx*kx+ky*ky))
    return np.fft.ifft2(np.fft.fft2(a)*transfer).real


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
    r_deg=np.hypot(xx-cx,yy-cy)*(FOV_DEG/WIDTH)
    ap=r_deg<=0.25
    ann=(r_deg>=0.75)&(r_deg<=1.50)
    bg=a[ann]
    med=float(np.median(bg))
    mad=float(np.median(np.abs(bg-med)))
    sig=1.4826*mad if mad>0 else float(np.std(bg))
    return float((np.mean(a[ap])-med)/sig) if sig>0 else None


def phase_markers(image: np.ndarray) -> dict[str,float]:
    a=finite_image(image)
    ny,nx=a.shape
    centered=a-float(np.median(a))
    F=np.fft.fftshift(np.fft.fft2(centered*np.outer(np.hanning(ny),np.hanning(nx))))
    amp=np.abs(F); power=amp*amp; phase=np.angle(F)
    yy,xx=np.indices(a.shape); cy,cx=(ny-1)/2.0,(nx-1)/2.0
    kx=(xx-cx)/nx; ky=(yy-cy)/ny; kr=np.hypot(kx,ky); th=np.arctan2(ky,kx)
    band=(kr>0.02)&(kr<0.20)
    aw=amp[band]; pw=power[band]
    coherence=float(abs(np.sum(aw*np.exp(1j*phase[band])))/np.sum(aw))
    hist,_=np.histogram(phase[band],bins=24,range=(-np.pi,np.pi),weights=aw)
    p=hist/hist.sum(); p=p[p>0]
    entropy=float(-np.sum(p*np.log(p))/np.log(24))
    out={"phase_coherence":coherence,"phase_entropy":entropy}
    for n in (1,2,3,6):
        out[f"A{n}"]=float(abs(np.sum(pw*np.exp(1j*n*th[band])))/np.sum(pw))
    return out


def phase36_sector_cell(image: np.ndarray) -> list[float]:
    a=finite_image(image)
    ny,nx=a.shape
    centered=a-float(np.median(a))
    F=np.fft.fftshift(np.fft.fft2(centered*np.outer(np.hanning(ny),np.hanning(nx))))
    amp=np.abs(F); ph=np.angle(F)
    yy,xx=np.indices(a.shape); cy,cx=(ny-1)/2.0,(nx-1)/2.0
    kx=(xx-cx)/nx; ky=(yy-cy)/ny; kr=np.hypot(kx,ky); th=np.arctan2(ky,kx)
    re=np.linspace(0.02,0.20,7); ae=np.linspace(-math.pi,math.pi,7)
    phases=[]
    for ri in range(6):
        for ai in range(6):
            m=(kr>=re[ri])&(kr<re[ri+1])&(th>=ae[ai])&(th<ae[ai+1])
            w=amp[m]; p=ph[m]
            if w.size==0 or float(np.sum(w))<=1e-18:
                phases.append(0.0)
            else:
                z=np.sum(w*np.exp(1j*p))
                phases.append(float(np.angle(z)%(2*math.pi)))
    if len(phases)!=36:
        raise RuntimeError("T36 sectorizer failed")
    return phases


def circular_delta(a: float,b: float) -> float:
    return math.atan2(math.sin(a-b),math.cos(a-b))


def moire_overlap(a: list[float],b: list[float]) -> float:
    z=sum(np.exp(1j*circular_delta(x,y)) for x,y in zip(a,b))/36.0
    return float(abs(z))


def mean_pairwise_overlap(cells: dict[str,list[float]]) -> float:
    names=sorted(cells)
    vals=[]
    for i,a in enumerate(names):
        for b in names[i+1:]:
            vals.append(moire_overlap(cells[a],cells[b]))
    return float(np.mean(vals))


def hologram_coherence(cells: dict[str,list[float]]) -> float:
    names=sorted(cells)
    amps=[]
    for j in range(36):
        z=sum(np.exp(1j*cells[name][j]) for name in names)
        amps.append(abs(z)/len(names))
    return float(np.mean(amps))


def fetch(ra_deg: float,dec_deg: float,source: dict[str,Any]) -> tuple[np.ndarray,str,str]:
    params={
        "hips":source["hips"],"width":str(WIDTH),"height":str(HEIGHT),
        "projection":"TAN","fov":str(FOV_DEG),"ra":str(ra_deg),"dec":str(dec_deg),
        "coordsys":"icrs","rotation_angle":"0","format":"fits"
    }
    url=HIPS2FITS+"?"+urlencode(params)
    r=requests.get(url,timeout=120,headers={"User-Agent":"GREMLIN-FERRET-PBH-NULL/0.2 read-only"})
    if r.status_code!=200 or not r.content.startswith(b"SIMPLE"):
        raise RuntimeError(f"bad FITS response {r.status_code}")
    with fits.open(io.BytesIO(r.content),memmap=False) as hdul:
        data=np.asarray(hdul[0].data,dtype=float)
        while data.ndim>2: data=data[0]
    return finite_image(data),hashlib.sha256(r.content).hexdigest(),url


def analyze_position(ra_deg: float,dec_deg: float,l_deg: float,b_deg: float) -> dict[str,Any]:
    raw={}; smooth={}; cells={}; hashes={}
    markers={}
    for name,src in SOURCES.items():
        data,sha,url=fetch(ra_deg,dec_deg,src)
        raw[name]=data
        sm=common_beam(data,src["native_fwhm_arcmin"])
        smooth[name]=sm
        q=WIDTH//4
        central=sm[q:3*q,q:3*q]
        cells[name]=phase36_sector_cell(central)
        markers[name]=phase_markers(central)
        hashes[name]=sha
    return {
        "ra_deg":ra_deg,"dec_deg":dec_deg,"l_deg":l_deg,"b_deg":b_deg,
        "center_z":{k:center_z(v) for k,v in smooth.items()},
        "marker":markers,
        "hashes":hashes,
        "smica_wmapV_image_corr":corr(smooth["planck_smica"],smooth["wmap_V"]),
        "smica_wmapW_image_corr":corr(smooth["planck_smica"],smooth["wmap_W"]),
        "wmapVW_image_corr":corr(smooth["wmap_V"],smooth["wmap_W"]),
        "gremlin_core_mean_pairwise_coherence":mean_pairwise_overlap(cells),
        "gremlin_hologram_coherence":hologram_coherence(cells),
    }


def circ_dist_deg(a: float,b: float) -> float:
    d=abs(a-b)%360.0
    return min(d,360.0-d)


def gal_to_icrs(l_deg: float,b_deg: float) -> tuple[float,float]:
    c=SkyCoord(l=l_deg*u.deg,b=b_deg*u.deg,frame="galactic").icrs
    return float(c.ra.deg),float(c.dec.deg)


def empirical_upper(target: float, controls: list[float]) -> float:
    n=sum(v>=target for v in controls)
    return (n+1)/(len(controls)+1)


def empirical_abs_upper(target: float, controls: list[float]) -> float:
    t=abs(target); n=sum(abs(v)>=t for v in controls)
    return (n+1)/(len(controls)+1)


def main() -> int:
    target=analyze_position(TARGET["ra_deg"],TARGET["dec_deg"],TARGET["gal_l_deg"],TARGET["gal_b_deg"])
    grid=[15.0*i for i in range(24)]
    control_l=[l for l in grid if circ_dist_deg(l,TARGET["gal_l_deg"])>=EXCLUSION_DEG]
    controls=[]
    for i,l in enumerate(control_l,1):
        ra,dec=gal_to_icrs(l,TARGET["gal_b_deg"])
        print(f"CONTROL {i}/{len(control_l)} l={l}",flush=True)
        controls.append(analyze_position(ra,dec,l,TARGET["gal_b_deg"]))

    coh=[x["gremlin_core_mean_pairwise_coherence"] for x in controls]
    holo=[x["gremlin_hologram_coherence"] for x in controls]
    zv=[x["center_z"]["planck_smica"] for x in controls]
    cv=[x["smica_wmapV_image_corr"] for x in controls]
    cw=[x["smica_wmapW_image_corr"] for x in controls]

    receipt={
        "schema":"PBH_CMB_MATCHED_LATITUDE_NULL_V0_2",
        "status":"ANALYZED",
        "target":target,
        "control_definition":{
            "gal_b_deg":TARGET["gal_b_deg"],
            "longitude_grid_deg":control_l,
            "exclusion_from_target_deg":EXCLUSION_DEG,
            "width":WIDTH,"height":HEIGHT,"fov_deg":FOV_DEG,
            "common_beam_fwhm_arcmin":COMMON_FWHM_ARCMIN,
        },
        "n_controls":len(controls),
        "controls_summary":[{
            "l_deg":x["l_deg"],"ra_deg":x["ra_deg"],"dec_deg":x["dec_deg"],
            "smica_center_z":x["center_z"]["planck_smica"],
            "smica_wmapV_image_corr":x["smica_wmapV_image_corr"],
            "smica_wmapW_image_corr":x["smica_wmapW_image_corr"],
            "gremlin_core_mean_pairwise_coherence":x["gremlin_core_mean_pairwise_coherence"],
            "gremlin_hologram_coherence":x["gremlin_hologram_coherence"],
        } for x in controls],
        "empirical":{
            "abs_smica_center_z":empirical_abs_upper(target["center_z"]["planck_smica"],zv),
            "gremlin_core_mean_pairwise_coherence_upper":empirical_upper(target["gremlin_core_mean_pairwise_coherence"],coh),
            "gremlin_hologram_coherence_upper":empirical_upper(target["gremlin_hologram_coherence"],holo),
            "smica_wmapV_image_corr_upper":empirical_upper(target["smica_wmapV_image_corr"],cv),
            "smica_wmapW_image_corr_upper":empirical_upper(target["smica_wmapW_image_corr"],cw),
        },
        "interpretation":"MATCHED_GALACTIC_LATITUDE_EMPIRICAL_NULL_NOT_FORMAL_SURVEY_PVALUE",
        "authority":"CANDIDATE_ONLY",
        "canon_allowed":False,
    }
    print("PBH_NULL_RECEIPT="+json.dumps(receipt,sort_keys=True,separators=(",",":")),flush=True)
    return 0

if __name__=="__main__":
    sys.exit(main())
