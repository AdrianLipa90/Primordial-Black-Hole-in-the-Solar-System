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

# Frozen before this validation run.
STATISTIC_ID="S_HIGH_MINUS_MID_T36_V0_1"
STATISTIC_FORMULA="mean(M353,M545,M857)-mean(M100,M143,M217)"
TAIL="UPPER"

SOURCES={
    "smica":{"hips":"CDS/P/PLANCK/R3/CMB","fwhm":5.0},
    "100":{"hips":"CDS/P/PLANCK/R3/HFI100","fwhm":9.66},
    "143":{"hips":"CDS/P/PLANCK/R3/HFI143","fwhm":7.27},
    "217":{"hips":"CDS/P/PLANCK/R3/HFI217","fwhm":5.01},
    "353":{"hips":"CDS/P/PLANCK/R3/HFI353","fwhm":4.86},
    "545":{"hips":"CDS/P/PLANCK/R3/HFI545","fwhm":4.84},
    "857":{"hips":"CDS/P/PLANCK/R3/HFI857","fwhm":4.63},
}

# Fresh grid: 7.5-degree offset from all earlier 15-degree control grids.
GRID_LONGITUDES=[7.5+15.0*i for i in range(24)]
EXCLUSION_DEG=30.0


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
    tf=np.exp(-2.0*(math.pi**2)*(sigma**2)*(kx*kx+ky*ky))
    return np.fft.ifft2(np.fft.fft2(a)*tf).real


def t36(a):
    a=finite(a); ny,nx=a.shape
    centered=a-float(np.median(a))
    F=np.fft.fftshift(np.fft.fft2(centered*np.outer(np.hanning(ny),np.hanning(nx))))
    amp=np.abs(F); ph=np.angle(F)
    yy,xx=np.indices(a.shape)
    cy,cx=(ny-1)/2,(nx-1)/2
    kx=(xx-cx)/nx; ky=(yy-cy)/ny
    kr=np.hypot(kx,ky); th=np.arctan2(ky,kx)
    re=np.linspace(0.02,0.20,7)
    ae=np.linspace(-math.pi,math.pi,7)
    out=[]
    for ri in range(6):
        for ai in range(6):
            mask=(kr>=re[ri])&(kr<re[ri+1])&(th>=ae[ai])&(th<ae[ai+1])
            w=amp[mask]; p=ph[mask]
            if w.size==0 or float(w.sum())<=1e-18:
                out.append(0.0)
            else:
                out.append(float(np.angle(np.sum(w*np.exp(1j*p)))%(2*math.pi)))
    if len(out)!=36: raise RuntimeError("T36 dimension error")
    return out


def moire(a,b):
    if len(a)!=36 or len(b)!=36: raise ValueError("T36 required")
    z=sum(np.exp(1j*math.atan2(math.sin(x-y),math.cos(x-y))) for x,y in zip(a,b))/36.0
    return float(abs(z))


def fetch_one(name,src,ra,dec):
    params={
        "hips":src["hips"],"width":str(WIDTH),"height":str(HEIGHT),
        "projection":"TAN","fov":str(FOV_DEG),
        "ra":f"{ra:.12f}","dec":f"{dec:.12f}",
        "coordsys":"icrs","rotation_angle":"0","format":"fits",
    }
    url=HIPS2FITS+"?"+urlencode(params)
    r=requests.get(url,timeout=30,headers={"User-Agent":"GREMLIN-FERRET-PBH-T36-BLIND/0.1 read-only"})
    if r.status_code!=200 or not r.content.startswith(b"SIMPLE"):
        raise RuntimeError(f"{name}: HTTP {r.status_code} bytes={len(r.content)}")
    with fits.open(io.BytesIO(r.content),memmap=False) as hdul:
        data=np.asarray(hdul[0].data,float)
        while data.ndim>2: data=data[0]
    if data.shape!=(HEIGHT,WIDTH):
        raise RuntimeError(f"{name}: unexpected shape {data.shape}")
    return name,finite(data)


def analyze(ra,dec):
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
    return {
        "M":M,
        "mid_mean":mid,
        "high_mean":high,
        "S_high_minus_mid":high-mid,
    }


def lonsep(a,b):
    return abs((a-b+180.0)%360.0-180.0)


def gal_to_icrs(l,b):
    c=SkyCoord(l=l*u.deg,b=b*u.deg,frame="galactic").icrs
    return float(c.ra.deg),float(c.dec.deg)


def main():
    prereg={
        "schema":"PBH_T36_BLIND_VALIDATION_PREREG_V0_1",
        "statistic_id":STATISTIC_ID,
        "formula":STATISTIC_FORMULA,
        "tail":TAIL,
        "common_beam_fwhm_arcmin":COMMON_FWHM_ARCMIN,
        "field_deg":FOV_DEG,
        "central_analysis_field_deg":4.0,
        "t36_layout":"RADIAL6_X_ANGULAR6_ROW_MAJOR",
        "fresh_control_grid_longitudes_deg":GRID_LONGITUDES,
        "same_galactic_latitude_deg":TARGET["gal_b_deg"],
        "target_exclusion_deg":EXCLUSION_DEG,
        "note":"Grid is offset by 7.5 deg from prior 15-deg grids.",
    }
    print("PBH_T36_BLIND_PREREG="+json.dumps(prereg,sort_keys=True,separators=(",",":")),flush=True)

    target=analyze(TARGET["ra_deg"],TARGET["dec_deg"])
    print("PBH_T36_BLIND_TARGET="+json.dumps(target,sort_keys=True,separators=(",",":")),flush=True)

    controls=[]
    for l in GRID_LONGITUDES:
        if lonsep(l,TARGET["gal_l_deg"])<EXCLUSION_DEG:
            continue
        ra,dec=gal_to_icrs(l,TARGET["gal_b_deg"])
        a=analyze(ra,dec)
        controls.append({"l_deg":l,"b_deg":TARGET["gal_b_deg"],"ra_deg":ra,"dec_deg":dec,**a})
        print(f"PBH_T36_BLIND_CONTROL {len(controls)} l={l:.1f} S={a['S_high_minus_mid']:.9f}",flush=True)

    vals=np.asarray([r["S_high_minus_mid"] for r in controls],float)
    rank_upper=float((1+np.sum(vals>=target["S_high_minus_mid"]))/(len(vals)+1))
    z=(target["S_high_minus_mid"]-float(vals.mean()))/float(vals.std(ddof=1)) if len(vals)>1 and vals.std(ddof=1)>0 else None

    receipt={
        "schema":"PBH_CMB_T36_FRESH_GRID_VALIDATION_V0_1",
        "status":"ANALYZED",
        "preregistration":prereg,
        "target":{**TARGET,**target},
        "n_controls":len(controls),
        "controls":controls,
        "control_distribution":{
            "min":float(vals.min()),"median":float(np.median(vals)),
            "mean":float(vals.mean()),"std_sample":float(vals.std(ddof=1)),
            "max":float(vals.max()),
        },
        "validation":{
            "upper_add_one_rank":rank_upper,
            "standardized_against_control_mean":z,
        },
        "interpretation":"FRESH_CONTROL_GRID_VALIDATION_OF_POSTHOC_STATISTIC_NOT_FORMAL_DISCOVERY_PVALUE",
        "authority":"CANDIDATE_ONLY",
        "canon_allowed":False,
    }
    compact={
        "schema":receipt["schema"],"status":"ANALYZED",
        "target_S":target["S_high_minus_mid"],
        "n_controls":len(controls),
        "control_distribution":receipt["control_distribution"],
        "upper_add_one_rank":rank_upper,
        "standardized_against_control_mean":z,
        "interpretation":receipt["interpretation"],
    }
    print("PBH_T36_BLIND_SUMMARY="+json.dumps(compact,sort_keys=True,separators=(",",":")),flush=True)
    print("PBH_T36_BLIND_RECEIPT="+json.dumps(receipt,sort_keys=True,separators=(",",":")),flush=True)
    return 0

if __name__=="__main__":
    raise SystemExit(main())
