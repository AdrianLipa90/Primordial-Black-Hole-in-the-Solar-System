#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import io
import json
import math
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Any
from urllib.parse import urlencode

import numpy as np
import requests
from astropy.io import fits

TARGET = {"ra_deg": 107.576646, "dec_deg": -25.907579}
WIDTH = 256
HEIGHT = 256
FOV_DEG = 8.0
PIX_ARCMIN = FOV_DEG * 60.0 / WIDTH
COMMON_FWHM_ARCMIN = 60.0
HIPS2FITS = "https://alasky.cds.unistra.fr/hips-image-services/hips2fits"

SOURCES = {
    "smica": {"hips":"CDS/P/PLANCK/R3/CMB","freq_ghz":None,"fwhm":5.0,"class":"CMB_COMPONENT"},
    "030": {"hips":"CDS/P/PLANCK/R3/LFI30","freq_ghz":30.0,"fwhm":32.29,"class":"PLANCK_FREQUENCY"},
    "044": {"hips":"CDS/P/PLANCK/R3/LFI44","freq_ghz":44.0,"fwhm":27.0,"class":"PLANCK_FREQUENCY"},
    "070": {"hips":"CDS/P/PLANCK/R3/LFI70","freq_ghz":70.0,"fwhm":13.21,"class":"PLANCK_FREQUENCY"},
    "100": {"hips":"CDS/P/PLANCK/R3/HFI100","freq_ghz":100.0,"fwhm":9.66,"class":"PLANCK_FREQUENCY"},
    "143": {"hips":"CDS/P/PLANCK/R3/HFI143","freq_ghz":143.0,"fwhm":7.27,"class":"PLANCK_FREQUENCY"},
    "217": {"hips":"CDS/P/PLANCK/R3/HFI217","freq_ghz":217.0,"fwhm":5.01,"class":"PLANCK_FREQUENCY"},
    "353": {"hips":"CDS/P/PLANCK/R3/HFI353","freq_ghz":353.0,"fwhm":4.86,"class":"PLANCK_FREQUENCY"},
    "545": {"hips":"CDS/P/PLANCK/R3/HFI545","freq_ghz":545.0,"fwhm":4.84,"class":"PLANCK_FREQUENCY"},
    "857": {"hips":"CDS/P/PLANCK/R3/HFI857","freq_ghz":857.0,"fwhm":4.63,"class":"PLANCK_FREQUENCY"},
}


def finite_image(image):
    image=np.asarray(image,float)
    m=np.isfinite(image)
    if not m.any(): raise RuntimeError("no finite pixels")
    med=float(np.nanmedian(image[m]))
    return np.where(m,image,med)


def common_beam(image,native_fwhm):
    a=finite_image(image)
    if native_fwhm>=COMMON_FWHM_ARCMIN: return a
    extra=math.sqrt(COMMON_FWHM_ARCMIN**2-native_fwhm**2)
    sigma=extra/(2.354820045*PIX_ARCMIN)
    ny,nx=a.shape
    fy=np.fft.fftfreq(ny); fx=np.fft.fftfreq(nx)
    ky,kx=np.meshgrid(fy,fx,indexing="ij")
    transfer=np.exp(-2.0*(math.pi**2)*(sigma**2)*(kx*kx+ky*ky))
    return np.fft.ifft2(np.fft.fft2(a)*transfer).real


def phase36(image):
    a=finite_image(image)
    centered=a-float(np.median(a))
    ny,nx=a.shape
    F=np.fft.fftshift(np.fft.fft2(centered*np.outer(np.hanning(ny),np.hanning(nx))))
    amp=np.abs(F); ph=np.angle(F)
    yy,xx=np.indices(a.shape); cy,cx=(ny-1)/2,(nx-1)/2
    kx=(xx-cx)/nx; ky=(yy-cy)/ny
    kr=np.hypot(kx,ky); th=np.arctan2(ky,kx)
    redges=np.linspace(0.02,0.20,7); aedges=np.linspace(-math.pi,math.pi,7)
    out=[]
    for ri in range(6):
        for ai in range(6):
            m=(kr>=redges[ri])&(kr<redges[ri+1])&(th>=aedges[ai])&(th<aedges[ai+1])
            w=amp[m]; p=ph[m]
            if w.size==0 or float(w.sum())<=1e-18: out.append(0.0)
            else:
                z=np.sum(w*np.exp(1j*p)); out.append(float(np.angle(z)%(2*math.pi)))
    if len(out)!=36: raise RuntimeError("phase36 != 36")
    return out


def moire(a,b):
    z=sum(np.exp(1j*math.atan2(math.sin(x-y),math.cos(x-y))) for x,y in zip(a,b))/36.0
    return float(abs(z))


def corr(a,b):
    x=(finite_image(a)-np.mean(a)).ravel(); y=(finite_image(b)-np.mean(b)).ravel()
    d=float(np.linalg.norm(x)*np.linalg.norm(y))
    return float(np.dot(x,y)/d) if d else None


def fetch(name,src):
    params={
        "hips":src["hips"],"width":str(WIDTH),"height":str(HEIGHT),
        "projection":"TAN","fov":str(FOV_DEG),
        "ra":str(TARGET["ra_deg"]),"dec":str(TARGET["dec_deg"]),
        "coordsys":"icrs","rotation_angle":"0","format":"fits"
    }
    url=HIPS2FITS+"?"+urlencode(params)
    r=requests.get(url,timeout=30,headers={"User-Agent":"GREMLIN-FERRET-PBH-SPECTRAL/0.1 read-only"})
    if r.status_code!=200 or not r.content.startswith(b"SIMPLE"):
        return name,None,{"error":f"HTTP {r.status_code} bytes={len(r.content)}","url":url}
    with fits.open(io.BytesIO(r.content),memmap=False) as hdul:
        data=np.asarray(hdul[0].data,float)
        while data.ndim>2: data=data[0]
        hdr=hdul[0].header
    if data.shape!=(HEIGHT,WIDTH):
        return name,None,{"error":f"shape {data.shape}","url":url}
    return name,finite_image(data),{
        "url":url,"sha256":hashlib.sha256(r.content).hexdigest(),"bytes":len(r.content),
        "bunit":hdr.get("BUNIT"),"history":[str(v) for v in hdr.get("HISTORY",[])]
    }


def main():
    raw={}; meta={}
    with ThreadPoolExecutor(max_workers=10) as ex:
        futs={ex.submit(fetch,n,s):n for n,s in SOURCES.items()}
        for fut in as_completed(futs):
            n,img,m=fut.result(); meta[n]=m
            print(f"SPECTRAL_FETCH {n} " + ("OK" if img is not None else "FAIL"),flush=True)
            if img is not None: raw[n]=img

    if "smica" not in raw: raise RuntimeError("SMICA unavailable")
    common={}
    cells={}
    q=WIDTH//4
    for n,img in raw.items():
        sm=common_beam(img,SOURCES[n]["fwhm"])
        common[n]=sm
        cells[n]=phase36(sm[q:3*q,q:3*q])

    smica=common["smica"]; smica_cell=cells["smica"]
    rows=[]
    for n,src in sorted(SOURCES.items(),key=lambda kv:(-1 if kv[1]["freq_ghz"] is None else kv[1]["freq_ghz"])):
        if n=="smica": continue
        if n not in common:
            rows.append({"id":n,"freq_ghz":src["freq_ghz"],"error":meta[n].get("error")})
            continue
        rows.append({
            "id":n,
            "freq_ghz":src["freq_ghz"],
            "smica_image_corr_common60":corr(smica,common[n]),
            "smica_t36_moire_common60":moire(smica_cell,cells[n]),
            "sha256":meta[n]["sha256"],
            "bunit":meta[n].get("bunit")
        })

    lows=[r for r in rows if r.get("freq_ghz") and r["freq_ghz"]<=70 and "error" not in r]
    mids=[r for r in rows if r.get("freq_ghz") and 100<=r["freq_ghz"]<=217 and "error" not in r]
    highs=[r for r in rows if r.get("freq_ghz") and r["freq_ghz"]>=353 and "error" not in r]
    def avg(group,key):
        return float(np.mean([r[key] for r in group])) if group else None

    receipt={
        "schema":"PBH_CMB_PLANCK_SPECTRAL_T36_GATE_V0_1",
        "status":"ANALYZED",
        "target":TARGET,
        "common_beam_fwhm_arcmin":COMMON_FWHM_ARCMIN,
        "rows":rows,
        "band_summary":{
            "low_30_70_t36_mean":avg(lows,"smica_t36_moire_common60"),
            "mid_100_217_t36_mean":avg(mids,"smica_t36_moire_common60"),
            "high_353_857_t36_mean":avg(highs,"smica_t36_moire_common60"),
            "low_30_70_image_corr_mean":avg(lows,"smica_image_corr_common60"),
            "mid_100_217_image_corr_mean":avg(mids,"smica_image_corr_common60"),
            "high_353_857_image_corr_mean":avg(highs,"smica_image_corr_common60"),
        },
        "interpretation":"FOREGROUND_SPECTRAL_DIAGNOSTIC_ONLY_NOT_PBH_DETECTION",
        "authority":"CANDIDATE_ONLY",
        "canon_allowed":False,
    }
    compact={
        "schema":receipt["schema"],"status":receipt["status"],
        "rows":rows,"band_summary":receipt["band_summary"],
        "interpretation":receipt["interpretation"]
    }
    print("PBH_SPECTRAL_T36_SUMMARY="+json.dumps(compact,sort_keys=True,separators=(",",":")),flush=True)
    print("PBH_SPECTRAL_T36_RECEIPT="+json.dumps(receipt,sort_keys=True,separators=(",",":")),flush=True)
    return 0


if __name__=="__main__":
    raise SystemExit(main())
