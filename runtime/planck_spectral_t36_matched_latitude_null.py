#!/usr/bin/env python3
from __future__ import annotations

import io, json, math
from concurrent.futures import ThreadPoolExecutor, as_completed
from urllib.parse import urlencode

import numpy as np
import requests
from astropy.coordinates import SkyCoord
from astropy.io import fits
import astropy.units as u

TARGET = {
    "ra_deg":107.576646,
    "dec_deg":-25.907579,
    "gal_l_deg":238.1709218507047,
    "gal_b_deg":-7.665782897020605,
}
WIDTH=256
HEIGHT=256
FOV_DEG=8.0
PIX_ARCMIN=FOV_DEG*60.0/WIDTH
COMMON_FWHM_ARCMIN=60.0
HIPS2FITS="https://alasky.cds.unistra.fr/hips-image-services/hips2fits"

SOURCES={
    "smica":{"hips":"CDS/P/PLANCK/R3/CMB","freq":None,"fwhm":5.0},
    "030":{"hips":"CDS/P/PLANCK/R3/LFI30","freq":30.0,"fwhm":32.29},
    "044":{"hips":"CDS/P/PLANCK/R3/LFI44","freq":44.0,"fwhm":27.0},
    "070":{"hips":"CDS/P/PLANCK/R3/LFI70","freq":70.0,"fwhm":13.21},
    "100":{"hips":"CDS/P/PLANCK/R3/HFI100","freq":100.0,"fwhm":9.66},
    "143":{"hips":"CDS/P/PLANCK/R3/HFI143","freq":143.0,"fwhm":7.27},
    "217":{"hips":"CDS/P/PLANCK/R3/HFI217","freq":217.0,"fwhm":5.01},
    "353":{"hips":"CDS/P/PLANCK/R3/HFI353","freq":353.0,"fwhm":4.86},
    "545":{"hips":"CDS/P/PLANCK/R3/HFI545","freq":545.0,"fwhm":4.84},
    "857":{"hips":"CDS/P/PLANCK/R3/HFI857","freq":857.0,"fwhm":4.63},
}
CONTROL_LONGITUDES=[0,15,30,45,60,75,90,105,120,135,150,165,180,195,270,285,300,315,330,345]


def finite(a):
    a=np.asarray(a,float); m=np.isfinite(a)
    if not m.any(): raise RuntimeError("no finite pixels")
    med=float(np.nanmedian(a[m]))
    return np.where(m,a,med)


def smooth(a,fwhm):
    a=finite(a)
    if fwhm>=COMMON_FWHM_ARCMIN: return a
    extra=math.sqrt(COMMON_FWHM_ARCMIN**2-fwhm**2)
    sigma=extra/(2.354820045*PIX_ARCMIN)
    ny,nx=a.shape; fy=np.fft.fftfreq(ny); fx=np.fft.fftfreq(nx)
    ky,kx=np.meshgrid(fy,fx,indexing="ij")
    tf=np.exp(-2*(math.pi**2)*(sigma**2)*(kx*kx+ky*ky))
    return np.fft.ifft2(np.fft.fft2(a)*tf).real


def t36(a):
    a=finite(a); ny,nx=a.shape
    F=np.fft.fftshift(np.fft.fft2((a-np.median(a))*np.outer(np.hanning(ny),np.hanning(nx))))
    amp=np.abs(F); ph=np.angle(F)
    yy,xx=np.indices(a.shape); cy,cx=(ny-1)/2,(nx-1)/2
    kx=(xx-cx)/nx; ky=(yy-cy)/ny; kr=np.hypot(kx,ky); th=np.arctan2(ky,kx)
    re=np.linspace(0.02,0.20,7); ae=np.linspace(-math.pi,math.pi,7)
    out=[]
    for ri in range(6):
        for ai in range(6):
            m=(kr>=re[ri])&(kr<re[ri+1])&(th>=ae[ai])&(th<ae[ai+1])
            w=amp[m]; p=ph[m]
            if not len(w) or float(w.sum())<=1e-18: out.append(0.0)
            else: out.append(float(np.angle(np.sum(w*np.exp(1j*p)))%(2*math.pi)))
    return out


def moire(a,b):
    z=sum(np.exp(1j*math.atan2(math.sin(x-y),math.cos(x-y))) for x,y in zip(a,b))/36.0
    return float(abs(z))


def fetch_one(name,src,ra,dec):
    params={"hips":src["hips"],"width":str(WIDTH),"height":str(HEIGHT),"projection":"TAN","fov":str(FOV_DEG),
            "ra":f"{ra:.12f}","dec":f"{dec:.12f}","coordsys":"icrs","rotation_angle":"0","format":"fits"}
    url=HIPS2FITS+"?"+urlencode(params)
    r=requests.get(url,timeout=30,headers={"User-Agent":"GREMLIN-FERRET-PBH-SPECTRAL-NULL/0.1 read-only"})
    if r.status_code!=200 or not r.content.startswith(b"SIMPLE"):
        raise RuntimeError(f"{name}: HTTP {r.status_code} bytes={len(r.content)}")
    with fits.open(io.BytesIO(r.content),memmap=False) as hdul:
        data=np.asarray(hdul[0].data,float)
        while data.ndim>2: data=data[0]
    if data.shape!=(HEIGHT,WIDTH): raise RuntimeError(f"{name}: shape {data.shape}")
    return name,finite(data)


def analyze(ra,dec):
    raw={}
    with ThreadPoolExecutor(max_workers=10) as ex:
        futs={ex.submit(fetch_one,n,s,ra,dec):n for n,s in SOURCES.items()}
        for fut in as_completed(futs):
            n,img=fut.result(); raw[n]=img
    q=WIDTH//4
    cells={}
    for n,src in SOURCES.items():
        sm=smooth(raw[n],src["fwhm"])
        cells[n]=t36(sm[q:3*q,q:3*q])
    s=cells["smica"]
    curve={n:moire(s,cells[n]) for n in SOURCES if n!="smica"}
    low=np.mean([curve["030"],curve["044"],curve["070"]])
    mid=np.mean([curve["100"],curve["143"],curve["217"]])
    high=np.mean([curve["353"],curve["545"],curve["857"]])
    peak70_margin=curve["070"]-max(curve["030"],curve["044"],curve["100"],curve["143"],curve["217"])
    return {
        "curve":curve,
        "t70":curve["070"],
        "low_mean":float(low),
        "mid_mean":float(mid),
        "high_mean":float(high),
        "peak70_margin":float(peak70_margin),
        "high_minus_mid":float(high-mid),
        "low_minus_mid":float(low-mid),
    }


def lonsep(a,b):
    return abs((a-b+180)%360-180)


def gal_to_icrs(l,b):
    c=SkyCoord(l=l*u.deg,b=b*u.deg,frame="galactic").icrs
    return float(c.ra.deg),float(c.dec.deg)


def upper_rank(v,controls):
    return float((1+sum(x>=v for x in controls))/(len(controls)+1))


def lower_rank(v,controls):
    return float((1+sum(x<=v for x in controls))/(len(controls)+1))


def main():
    target=analyze(TARGET["ra_deg"],TARGET["dec_deg"])
    print("SPECTRAL_NULL_TARGET="+json.dumps(target,sort_keys=True,separators=(",",":")),flush=True)
    controls=[]
    for l in CONTROL_LONGITUDES:
        if lonsep(l,TARGET["gal_l_deg"])<30: continue
        ra,dec=gal_to_icrs(l,TARGET["gal_b_deg"])
        a=analyze(ra,dec)
        row={"l_deg":l,"b_deg":TARGET["gal_b_deg"],"ra_deg":ra,"dec_deg":dec,**a}
        controls.append(row)
        print(f"SPECTRAL_NULL_CONTROL {len(controls)} l={l} t70={a['t70']:.6f} margin={a['peak70_margin']:.6f}",flush=True)

    fields=["t70","low_mean","mid_mean","high_mean","peak70_margin","high_minus_mid","low_minus_mid"]
    empirical={}
    distributions={}
    for field in fields:
        vals=[r[field] for r in controls]
        empirical[field+"_upper"]=upper_rank(target[field],vals)
        empirical[field+"_lower"]=lower_rank(target[field],vals)
        distributions[field]={
            "min":float(np.min(vals)),"median":float(np.median(vals)),
            "mean":float(np.mean(vals)),"max":float(np.max(vals))
        }

    receipt={
        "schema":"PBH_CMB_PLANCK_SPECTRAL_T36_MATCHED_LATITUDE_NULL_V0_1",
        "status":"ANALYZED",
        "target":{**TARGET,**target},
        "n_controls":len(controls),
        "controls":controls,
        "control_definition":{
            "same_galactic_latitude_deg":TARGET["gal_b_deg"],
            "target_exclusion_deg":30.0,
            "longitudes_deg":CONTROL_LONGITUDES,
            "common_beam_fwhm_arcmin":COMMON_FWHM_ARCMIN,
            "frequencies_ghz":[30,44,70,100,143,217,353,545,857],
        },
        "empirical_descriptive_ranks":empirical,
        "control_distributions":distributions,
        "interpretation":"POSTHOC_SPECTRAL_FOREGROUND_DIAGNOSTIC_ONLY_NOT_FORMAL_PVALUE_OR_PBH_DETECTION",
        "authority":"CANDIDATE_ONLY","canon_allowed":False
    }
    compact={
        "schema":receipt["schema"],"status":"ANALYZED","n_controls":len(controls),
        "target":{k:target[k] for k in fields},
        "empirical_descriptive_ranks":empirical,
        "control_distributions":distributions,
        "interpretation":receipt["interpretation"]
    }
    print("PBH_SPECTRAL_T36_NULL_SUMMARY="+json.dumps(compact,sort_keys=True,separators=(",",":")),flush=True)
    print("PBH_SPECTRAL_T36_NULL_RECEIPT="+json.dumps(receipt,sort_keys=True,separators=(",",":")),flush=True)
    return 0

if __name__=="__main__":
    raise SystemExit(main())
