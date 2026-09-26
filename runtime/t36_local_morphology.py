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

TARGET = {"gal_l_deg":238.1709218507047,"gal_b_deg":-7.665782897020605}
OFFSETS=[-6.0,0.0,6.0]
WIDTH=256
HEIGHT=256
FOV_DEG=8.0
PIX_ARCMIN=FOV_DEG*60.0/WIDTH
COMMON_FWHM_ARCMIN=60.0
HIPS2FITS="https://alasky.cds.unistra.fr/hips-image-services/hips2fits"

SOURCES={
 "smica":{"hips":"CDS/P/PLANCK/R3/CMB","fwhm":5.0},
 "100":{"hips":"CDS/P/PLANCK/R3/HFI100","fwhm":9.66},
 "143":{"hips":"CDS/P/PLANCK/R3/HFI143","fwhm":7.27},
 "217":{"hips":"CDS/P/PLANCK/R3/HFI217","fwhm":5.01},
 "353":{"hips":"CDS/P/PLANCK/R3/HFI353","fwhm":4.86},
 "545":{"hips":"CDS/P/PLANCK/R3/HFI545","fwhm":4.84},
 "857":{"hips":"CDS/P/PLANCK/R3/HFI857","fwhm":4.63},
}

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
    ny,nx=a.shape
    fy=np.fft.fftfreq(ny); fx=np.fft.fftfreq(nx)
    ky,kx=np.meshgrid(fy,fx,indexing="ij")
    tf=np.exp(-2*(math.pi**2)*(sigma**2)*(kx*kx+ky*ky))
    return np.fft.ifft2(np.fft.fft2(a)*tf).real

def t36(a):
    a=finite(a); ny,nx=a.shape
    F=np.fft.fftshift(np.fft.fft2((a-np.median(a))*np.outer(np.hanning(ny),np.hanning(nx))))
    amp=np.abs(F); ph=np.angle(F)
    yy,xx=np.indices(a.shape); cy,cx=(ny-1)/2,(nx-1)/2
    kx=(xx-cx)/nx; ky=(yy-cy)/ny
    kr=np.hypot(kx,ky); th=np.arctan2(ky,kx)
    re=np.linspace(0.02,0.20,7); ae=np.linspace(-math.pi,math.pi,7)
    out=[]
    for ri in range(6):
        for ai in range(6):
            m=(kr>=re[ri])&(kr<re[ri+1])&(th>=ae[ai])&(th<ae[ai+1])
            w=amp[m]; p=ph[m]
            if w.size==0 or float(w.sum())<=1e-18: out.append(0.0)
            else: out.append(float(np.angle(np.sum(w*np.exp(1j*p)))%(2*math.pi)))
    if len(out)!=36: raise RuntimeError("T36 dimension error")
    return out

def moire(a,b):
    z=sum(np.exp(1j*math.atan2(math.sin(x-y),math.cos(x-y))) for x,y in zip(a,b))/36.0
    return float(abs(z))

def gal_to_icrs(l,b):
    c=SkyCoord(l=l*u.deg,b=b*u.deg,frame="galactic").icrs
    return float(c.ra.deg),float(c.dec.deg)

def fetch_one(name,src,ra,dec):
    params={"hips":src["hips"],"width":str(WIDTH),"height":str(HEIGHT),
            "projection":"TAN","fov":str(FOV_DEG),"ra":f"{ra:.12f}",
            "dec":f"{dec:.12f}","coordsys":"icrs","rotation_angle":"0","format":"fits"}
    url=HIPS2FITS+"?"+urlencode(params)
    r=requests.get(url,timeout=30,headers={"User-Agent":"GREMLIN-FERRET-PBH-LOCAL-MORPH/0.1 read-only"})
    if r.status_code!=200 or not r.content.startswith(b"SIMPLE"):
        raise RuntimeError(f"{name}: HTTP {r.status_code} bytes={len(r.content)}")
    with fits.open(io.BytesIO(r.content),memmap=False) as hdul:
        data=np.asarray(hdul[0].data,float)
        while data.ndim>2: data=data[0]
    if data.shape!=(HEIGHT,WIDTH): raise RuntimeError(f"{name}: shape {data.shape}")
    return name,finite(data)

def analyze(l,b):
    ra,dec=gal_to_icrs(l,b)
    raw={}
    with ThreadPoolExecutor(max_workers=7) as ex:
        futs={ex.submit(fetch_one,n,s,ra,dec):n for n,s in SOURCES.items()}
        for fut in as_completed(futs):
            n,img=fut.result(); raw[n]=img
    q=WIDTH//4
    cells={}
    for n,src in SOURCES.items():
        sm=smooth(raw[n],src["fwhm"])
        cells[n]=t36(sm[q:3*q,q:3*q])
    ref=cells["smica"]
    M={n:moire(ref,cells[n]) for n in ("100","143","217","353","545","857")}
    mid=float(np.mean([M["100"],M["143"],M["217"]]))
    high=float(np.mean([M["353"],M["545"],M["857"]]))
    return {"l_deg":l,"b_deg":b,"ra_deg":ra,"dec_deg":dec,"mid_mean":mid,"high_mean":high,"S":high-mid,"M":M}

def main():
    prereg={
      "schema":"PBH_T36_LOCAL_MORPH_PREREG_V0_1",
      "statistic":"S=mean(M353,M545,M857)-mean(M100,M143,M217)",
      "grid_offsets_deg":OFFSETS,
      "coordinate_frame":"galactic",
      "common_beam_fwhm_arcmin":COMMON_FWHM_ARCMIN,
      "field_deg":FOV_DEG,
      "central_analysis_field_deg":4.0,
      "interpretation":"LOCALIZATION_DIAGNOSTIC_ONLY",
    }
    print("PBH_T36_LOCAL_PREREG="+json.dumps(prereg,sort_keys=True,separators=(",",":")),flush=True)
    rows=[]
    for db in OFFSETS:
        for dl in OFFSETS:
            l=(TARGET["gal_l_deg"]+dl)%360.0
            b=TARGET["gal_b_deg"]+db
            row=analyze(l,b)
            row["dl_deg"]=dl; row["db_deg"]=db
            rows.append(row)
            print(f"PBH_T36_LOCAL dl={dl:+.1f} db={db:+.1f} S={row['S']:.9f}",flush=True)

    center=next(r for r in rows if r["dl_deg"]==0 and r["db_deg"]==0)
    vals=np.asarray([r["S"] for r in rows],float)
    rank_desc=int(1+np.sum(vals>center["S"]))
    receipt={
      "schema":"PBH_CMB_T36_LOCAL_MORPHOLOGY_V0_1",
      "status":"ANALYZED",
      "preregistration":prereg,
      "target":TARGET,
      "rows":rows,
      "center_S":center["S"],
      "center_rank_descending_of_9":rank_desc,
      "grid_distribution":{
        "min":float(vals.min()),"median":float(np.median(vals)),
        "mean":float(vals.mean()),"max":float(vals.max())
      },
      "interpretation":"CENTER_LOCALIZATION_DIAGNOSTIC_ONLY_NOT_PBH_DETECTION",
      "authority":"CANDIDATE_ONLY","canon_allowed":False
    }
    compact={k:receipt[k] for k in ("schema","status","center_S","center_rank_descending_of_9","grid_distribution","interpretation")}
    print("PBH_T36_LOCAL_SUMMARY="+json.dumps(compact,sort_keys=True,separators=(",",":")),flush=True)
    print("PBH_T36_LOCAL_RECEIPT="+json.dumps(receipt,sort_keys=True,separators=(",",":")),flush=True)
    return 0

if __name__=="__main__":
    raise SystemExit(main())
