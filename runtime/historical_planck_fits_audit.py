#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import io
import json
import math
from itertools import combinations
from urllib.parse import quote

import numpy as np
import requests
from astropy.io import fits
from astropy.wcs import WCS

TARGET_RA_DEG = 107.57664605896267
TARGET_DEC_DEG = -25.907578906023748
SOURCE_REPO = "AdrianLipa90/Metatime"
SOURCE_COMMIT = "aa9fe46d2ca0a93e804e78ea843bbc296aeb591d"

FILES = {
    "planck_100_R3": (
        "107_-26_Equ_J2000_planck_100Rel3.fits",
        "1bef9697e314ad212d8aa2d1eec606968e189ead61a2d07e3d9eebec691409c5",
    ),
    "planck_143_R1": (
        "107_-26_Equ_J2000_planck_143Rel1.fits",
        "893e21091e6c9b770c9acc17a88ae874bc57ccef8a0e3dd35c5bce5552de1f8b",
    ),
    "planck_143_R2": (
        "107_-26_Equ_J2000_planck_143Rel2.fits",
        "63ebab928422961985d69b537744785d8321e5f023dfe523bb6586be524f1748",
    ),
    "planck_143_R3": (
        "107_-26_Equ_J2000_planck_143Rel3.fits",
        "61770ad10ac6eb20b85a6587da21149b5e46d69b33273ff9f5eea93ed09d44a1",
    ),
    "planck_217_R3": (
        "107_-26_Equ_J2000_planck_217Rel3.fits",
        "ba03fb5d2aaedfc56de6079529cd6bdc22ef875dd5019d6d2abcb7c5c0e7cd1a",
    ),
}


def finite_image(a):
    a=np.asarray(a,dtype=float)
    m=np.isfinite(a)
    if not m.any():
        raise RuntimeError("no finite pixels")
    med=float(np.nanmedian(a[m]))
    return np.where(m,a,med)


def corr(a,b):
    if a.shape != b.shape:
        return None
    x=(finite_image(a)-np.mean(a)).ravel()
    y=(finite_image(b)-np.mean(b)).ravel()
    d=float(np.linalg.norm(x)*np.linalg.norm(y))
    return float(np.dot(x,y)/d) if d else None


def phase_markers(image):
    a=finite_image(image)
    ny,nx=a.shape
    centered=a-float(np.median(a))
    F=np.fft.fftshift(np.fft.fft2(centered*np.outer(np.hanning(ny),np.hanning(nx))))
    amp=np.abs(F); power=amp*amp; ph=np.angle(F)
    yy,xx=np.indices(a.shape); cy,cx=(ny-1)/2,(nx-1)/2
    kx=(xx-cx)/nx; ky=(yy-cy)/ny; kr=np.hypot(kx,ky); th=np.arctan2(ky,kx)
    band=(kr>0.02)&(kr<0.20)
    aw=amp[band]; pw=power[band]
    if not band.any() or float(np.sum(aw))<=0:
        return {"error":"EMPTY_PHASE_BAND"}
    coh=float(abs(np.sum(aw*np.exp(1j*ph[band])))/np.sum(aw))
    hist,_=np.histogram(ph[band],bins=24,range=(-np.pi,np.pi),weights=aw)
    p=hist/hist.sum(); p=p[p>0]
    ent=float(-np.sum(p*np.log(p))/np.log(24))
    out={"phase_coherence":coh,"phase_entropy":ent}
    for n in (1,2,3,6):
        out[f"A{n}"]=float(abs(np.sum(pw*np.exp(1j*n*th[band])))/np.sum(pw))
    return out


def phase36(image):
    a=finite_image(image)
    ny,nx=a.shape
    centered=a-float(np.median(a))
    F=np.fft.fftshift(np.fft.fft2(centered*np.outer(np.hanning(ny),np.hanning(nx))))
    amp=np.abs(F); ph=np.angle(F)
    yy,xx=np.indices(a.shape); cy,cx=(ny-1)/2,(nx-1)/2
    kx=(xx-cx)/nx; ky=(yy-cy)/ny; kr=np.hypot(kx,ky); th=np.arctan2(ky,kx)
    re=np.linspace(0.02,0.20,7); ae=np.linspace(-math.pi,math.pi,7)
    vec=[]
    for ri in range(6):
        for ai in range(6):
            mask=(kr>=re[ri])&(kr<re[ri+1])&(th>=ae[ai])&(th<ae[ai+1])
            w=amp[mask]; p=ph[mask]
            if w.size==0 or float(np.sum(w))<=1e-18:
                vec.append(0.0)
            else:
                vec.append(float(np.angle(np.sum(w*np.exp(1j*p)))%(2*math.pi)))
    if len(vec)!=36:
        raise RuntimeError("T36 shape failure")
    return vec


def moire(a,b):
    z=sum(np.exp(1j*math.atan2(math.sin(x-y),math.cos(x-y))) for x,y in zip(a,b))/36.0
    return float(abs(z))


def target_patch(data, header, size_deg=4.0):
    w=WCS(header)
    px,py=w.all_world2pix(TARGET_RA_DEG,TARGET_DEC_DEG,0)
    ix=int(round(float(px))); iy=int(round(float(py)))
    cd1=abs(float(header.get("CDELT1",0.0)))
    cd2=abs(float(header.get("CDELT2",0.0)))
    pixscale=max(cd1,cd2)
    if pixscale<=0:
        raise RuntimeError("missing CDELT pixel scale")
    half=max(16,int(round((size_deg/2.0)/pixscale)))
    y0=max(0,iy-half); y1=min(data.shape[0],iy+half)
    x0=max(0,ix-half); x1=min(data.shape[1],ix+half)
    patch=data[y0:y1,x0:x1]
    if min(patch.shape)<32:
        raise RuntimeError(f"target patch too small {patch.shape}")
    return patch,(float(px),float(py)),(x0,x1,y0,y1),pixscale


def center_stats(data, header):
    w=WCS(header)
    px,py=w.all_world2pix(TARGET_RA_DEG,TARGET_DEC_DEG,0)
    yy,xx=np.indices(data.shape)
    pixscale=max(abs(float(header["CDELT1"])),abs(float(header["CDELT2"])))
    rdeg=np.hypot(xx-px,yy-py)*pixscale
    ap=rdeg<=0.25
    ann=(rdeg>=0.75)&(rdeg<=1.5)
    bg=data[ann]
    med=float(np.median(bg)); mad=float(np.median(np.abs(bg-med)))
    sig=1.4826*mad if mad>0 else float(np.std(bg))
    cm=float(np.mean(data[ap]))
    return {"center_mean":cm,"background_median":med,"robust_sigma":sig,"center_robust_z":float((cm-med)/sig) if sig>0 else None}


def main():
    receipt={
      "schema":"PBH_HISTORICAL_PLANCK_FITS_AUDIT_V0_1",
      "status":"PARTIAL",
      "source":{"repository":SOURCE_REPO,"commit":SOURCE_COMMIT},
      "target":{"ra_deg":TARGET_RA_DEG,"dec_deg":TARGET_DEC_DEG},
      "maps":{},"pairwise_image_correlations":{},"pairwise_t36_moire":{},
      "interpretation":"SOURCE_RECONSTRUCTION_ONLY_NOT_PBH_DETECTION",
      "authority":"CANDIDATE_ONLY","canon_allowed":False
    }
    images={}; cells={}
    for key,(name,expected_sha) in FILES.items():
        url=f"https://raw.githubusercontent.com/{SOURCE_REPO}/{SOURCE_COMMIT}/data/{quote(name)}"
        entry={"file":name,"url":url,"expected_sha256":expected_sha}
        try:
            r=requests.get(url,timeout=120,headers={"User-Agent":"FERRET-PBH-historical-FITS-audit/0.1"})
            entry["http_status"]=r.status_code; entry["bytes"]=len(r.content)
            actual=hashlib.sha256(r.content).hexdigest()
            entry["sha256"]=actual; entry["sha256_match"]=(actual==expected_sha)
            if r.status_code!=200 or not r.content.startswith(b"SIMPLE"):
                raise RuntimeError(f"invalid FITS response HTTP {r.status_code}")
            with fits.open(io.BytesIO(r.content),memmap=False) as hdul:
                h=hdul[0].header
                data=np.asarray(hdul[0].data,dtype=float)
                while data.ndim>2: data=data[0]
            data=finite_image(data)
            entry["shape"]=[int(data.shape[0]),int(data.shape[1])]
            entry["header"]={k:h.get(k) for k in (
              "SURVEY","TELESCOP","INSTRUME","BUNIT","CTYPE1","CTYPE2",
              "CRVAL1","CRVAL2","CRPIX1","CRPIX2","CDELT1","CDELT2",
              "CD1_1","CD1_2","CD2_1","CD2_2","RADESYS","RADECSYS","EQUINOX",
              "NAXIS1","NAXIS2"
            ) if k in h}
            patch,pix,bounds,pixscale=target_patch(data,h)
            entry.update({
              "target_pixel":[pix[0],pix[1]],"target_patch_bounds":list(bounds),"pixel_scale_deg":pixscale,
              "center_stats":center_stats(data,h),
              "phase_markers_target_4deg":phase_markers(patch),
              "t36_target_4deg":phase36(patch)
            })
            images[key]=data; cells[key]=entry["t36_target_4deg"]
        except Exception as exc:
            entry["error"]=f"{type(exc).__name__}: {exc}"
        receipt["maps"][key]=entry

    for a,b in combinations(sorted(images),2):
        receipt["pairwise_image_correlations"][f"{a}__{b}"]=corr(images[a],images[b])
        receipt["pairwise_t36_moire"][f"{a}__{b}"]=moire(cells[a],cells[b])

    rels=["planck_143_R1","planck_143_R2","planck_143_R3"]
    if all(x in images for x in rels):
        receipt["release_143_summary"]={
          "image_correlations":{
             "R1_R2":corr(images["planck_143_R1"],images["planck_143_R2"]),
             "R1_R3":corr(images["planck_143_R1"],images["planck_143_R3"]),
             "R2_R3":corr(images["planck_143_R2"],images["planck_143_R3"])
          },
          "t36_moire":{
             "R1_R2":moire(cells["planck_143_R1"],cells["planck_143_R2"]),
             "R1_R3":moire(cells["planck_143_R1"],cells["planck_143_R3"]),
             "R2_R3":moire(cells["planck_143_R2"],cells["planck_143_R3"])
          }
        }
    ok=[k for k,v in receipt["maps"].items() if "error" not in v and v.get("sha256_match")]
    receipt["verified_maps"]=ok
    receipt["status"]="ANALYZED" if len(ok)==len(FILES) else "PARTIAL"
    print("PBH_HIST_FITS_RECEIPT="+json.dumps(receipt,sort_keys=True,separators=(",",":")),flush=True)

if __name__=="__main__":
    main()
