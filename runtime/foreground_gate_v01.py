#!/usr/bin/env python3
from __future__ import annotations

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

TARGET={"ra_deg":107.576646,"dec_deg":-25.907579,"gal_l_deg":238.1709218507047,"gal_b_deg":-7.665782897020605}
WIDTH=128
HEIGHT=128
FOV_DEG=6.0
PIX_ARCMIN=FOV_DEG*60.0/WIDTH
COMMON_FWHM_ARCMIN=40.0
EXCLUSION_DEG=30.0
HIPS2FITS="https://alasky.cds.unistra.fr/hips-image-services/hips2fits"

SOURCES={
    "planck_smica":{"hips":"CDS/P/PLANCK/R3/CMB","native_fwhm_arcmin":5.0,"role":"CMB_COMPONENT"},
    "planck_353":{"hips":"CDS/P/PLANCK/R3/HFI353","native_fwhm_arcmin":4.94,"role":"DUST_SENSITIVE_CONTROL"},
    "planck_30":{"hips":"CDS/P/PLANCK/R3/LFI30","native_fwhm_arcmin":32.3,"role":"LOW_FREQUENCY_FOREGROUND_CONTROL"},
}


def finite_image(a):
    a=np.asarray(a,float)
    m=np.isfinite(a)
    if not m.any(): raise RuntimeError("no finite pixels")
    med=float(np.nanmedian(a[m]))
    return np.where(m,a,med)


def common_beam(image,native):
    a=finite_image(image)
    if native>=COMMON_FWHM_ARCMIN: return a
    extra=math.sqrt(COMMON_FWHM_ARCMIN**2-native**2)
    sigma=extra/(2.354820045*PIX_ARCMIN)
    ny,nx=a.shape
    fy=np.fft.fftfreq(ny); fx=np.fft.fftfreq(nx)
    ky,kx=np.meshgrid(fy,fx,indexing="ij")
    transfer=np.exp(-2.0*(math.pi**2)*(sigma**2)*(kx*kx+ky*ky))
    return np.fft.ifft2(np.fft.fft2(a)*transfer).real


def corr(a,b):
    x=(finite_image(a)-np.mean(a)).ravel()
    y=(finite_image(b)-np.mean(b)).ravel()
    d=float(np.linalg.norm(x)*np.linalg.norm(y))
    return float(np.dot(x,y)/d) if d>0 else None


def phase36(image):
    a=finite_image(image)
    ny,nx=a.shape
    centered=a-float(np.median(a))
    F=np.fft.fftshift(np.fft.fft2(centered*np.outer(np.hanning(ny),np.hanning(nx))))
    amp=np.abs(F); ph=np.angle(F)
    yy,xx=np.indices(a.shape); cy,cx=(ny-1)/2,(nx-1)/2
    kx=(xx-cx)/nx; ky=(yy-cy)/ny; kr=np.hypot(kx,ky); th=np.arctan2(ky,kx)
    re=np.linspace(0.02,0.20,7); ae=np.linspace(-math.pi,math.pi,7)
    out=[]
    for ri in range(6):
        for ai in range(6):
            m=(kr>=re[ri])&(kr<re[ri+1])&(th>=ae[ai])&(th<ae[ai+1])
            w=amp[m]; p=ph[m]
            if w.size==0 or float(np.sum(w))<=1e-18: out.append(0.0)
            else: out.append(float(np.angle(np.sum(w*np.exp(1j*p)))%(2*math.pi)))
    if len(out)!=36: raise RuntimeError("phase36 failed")
    return out


def delta(a,b): return math.atan2(math.sin(a-b),math.cos(a-b))


def moire(a,b):
    return float(abs(sum(np.exp(1j*delta(x,y)) for x,y in zip(a,b))/36.0))


def fetch(ra,dec,src):
    params={"hips":src["hips"],"width":str(WIDTH),"height":str(HEIGHT),"projection":"TAN","fov":str(FOV_DEG),"ra":str(ra),"dec":str(dec),"coordsys":"icrs","rotation_angle":"0","format":"fits"}
    url=HIPS2FITS+"?"+urlencode(params)
    r=requests.get(url,timeout=120,headers={"User-Agent":"GREMLIN-FERRET-PBH-FOREGROUND/0.1 read-only"})
    if r.status_code!=200 or not r.content.startswith(b"SIMPLE"):
        raise RuntimeError(f"bad FITS {r.status_code} {src['hips']}")
    with fits.open(io.BytesIO(r.content),memmap=False) as hdul:
        a=np.asarray(hdul[0].data,float)
        while a.ndim>2:a=a[0]
    return finite_image(a)


def analyze(ra,dec,l,b):
    maps={}
    cells={}
    for name,src in SOURCES.items():
        sm=common_beam(fetch(ra,dec,src),src["native_fwhm_arcmin"])
        maps[name]=sm
        q=WIDTH//4
        cells[name]=phase36(sm[q:3*q,q:3*q])
    return {
        "ra_deg":ra,"dec_deg":dec,"l_deg":l,"b_deg":b,
        "smica_353_image_corr":corr(maps["planck_smica"],maps["planck_353"]),
        "smica_30_image_corr":corr(maps["planck_smica"],maps["planck_30"]),
        "smica_353_t36_moire":moire(cells["planck_smica"],cells["planck_353"]),
        "smica_30_t36_moire":moire(cells["planck_smica"],cells["planck_30"]),
    }


def circ_dist(a,b):
    d=abs(a-b)%360.0
    return min(d,360.0-d)


def gal_to_icrs(l,b):
    c=SkyCoord(l=l*u.deg,b=b*u.deg,frame="galactic").icrs
    return float(c.ra.deg),float(c.dec.deg)


def p_upper(target,vals):
    return (sum(v>=target for v in vals)+1)/(len(vals)+1)


def p_abs_upper(target,vals):
    return (sum(abs(v)>=abs(target) for v in vals)+1)/(len(vals)+1)


def main():
    target=analyze(TARGET["ra_deg"],TARGET["dec_deg"],TARGET["gal_l_deg"],TARGET["gal_b_deg"])
    grid=[15.0*i for i in range(24)]
    control_l=[l for l in grid if circ_dist(l,TARGET["gal_l_deg"])>=EXCLUSION_DEG]
    controls=[]
    for i,l in enumerate(control_l,1):
        ra,dec=gal_to_icrs(l,TARGET["gal_b_deg"])
        print(f"FOREGROUND CONTROL {i}/{len(control_l)} l={l}",flush=True)
        controls.append(analyze(ra,dec,l,TARGET["gal_b_deg"]))
    metrics=["smica_353_image_corr","smica_30_image_corr","smica_353_t36_moire","smica_30_t36_moire"]
    empirical={}
    for k in metrics:
        vals=[x[k] for x in controls]
        if "image_corr" in k:
            empirical[k+"_abs_upper"]=p_abs_upper(target[k],vals)
        else:
            empirical[k+"_upper"]=p_upper(target[k],vals)
    receipt={
        "schema":"PBH_CMB_FOREGROUND_GATE_V0_1",
        "status":"ANALYZED",
        "target":target,
        "n_controls":len(controls),
        "controls":controls,
        "empirical":empirical,
        "control_definition":{"same_galactic_latitude":TARGET["gal_b_deg"],"longitude_grid_deg":control_l,"target_exclusion_deg":EXCLUSION_DEG,"common_beam_fwhm_arcmin":COMMON_FWHM_ARCMIN},
        "interpretation":"FOREGROUND_DIAGNOSTIC_ONLY_NOT_COMPONENT_SEPARATION_PROOF",
        "authority":"CANDIDATE_ONLY",
        "canon_allowed":False,
    }
    print("PBH_FOREGROUND_RECEIPT="+json.dumps(receipt,sort_keys=True,separators=(",",":")),flush=True)
    return 0

if __name__=="__main__": sys.exit(main())
