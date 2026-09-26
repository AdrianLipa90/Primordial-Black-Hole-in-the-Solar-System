#!/usr/bin/env python3
from __future__ import annotations
import argparse, hashlib, json
from pathlib import Path
import numpy as np
from astropy.io import fits

def sha256_file(path):
    h=hashlib.sha256()
    with open(path,"rb") as f:
        for chunk in iter(lambda:f.read(1024*1024),b""): h.update(chunk)
    return h.hexdigest()

def load_2d(path):
    with fits.open(path,memmap=False) as hdul:
        for hdu in hdul:
            if hdu.data is None: continue
            a=np.asarray(hdu.data)
            while a.ndim>2: a=a[0]
            if a.ndim==2: return np.asarray(a,float),dict(hdu.header)
    raise RuntimeError(f"no 2D image HDU: {path}")

def clean(a):
    m=np.isfinite(a)
    if not m.any(): raise RuntimeError("no finite pixels")
    med=float(np.nanmedian(a[m]))
    return np.where(m,a,med)

def markers(image):
    a=clean(image)
    mean=float(np.mean(a)); std=float(np.std(a))
    centered=a-float(np.median(a))
    ny,nx=centered.shape
    tapered=centered*np.outer(np.hanning(ny),np.hanning(nx))
    F=np.fft.fftshift(np.fft.fft2(tapered))
    amp=np.abs(F); power=amp*amp; phase=np.angle(F)
    y,x=np.indices(a.shape); cy,cx=(ny-1)/2,(nx-1)/2
    kx=(x-cx)/nx; ky=(y-cy)/ny; kr=np.hypot(kx,ky); th=np.arctan2(ky,kx)
    band=(kr>0.02)&(kr<0.20)
    aw=amp[band]; pw=power[band]
    coherence=float(abs(np.sum(aw*np.exp(1j*phase[band])))/np.sum(aw))
    hist,_=np.histogram(phase[band],bins=24,range=(-np.pi,np.pi),weights=aw)
    p=hist/hist.sum(); p=p[p>0]
    entropy=float(-np.sum(p*np.log(p))/np.log(24))
    out={"shape":[ny,nx],"mean":mean,"std":std,"phase_coherence":coherence,"phase_entropy":entropy}
    for n in (1,2,3,6):
        out[f"A{n}"]=float(abs(np.sum(pw*np.exp(1j*n*th[band])))/np.sum(pw))
    dom=(kr>0.01)&(kr<0.25); coords=np.argwhere(dom); yy,xx=coords[int(np.argmax(power[dom]))]
    out.update({"dominant_phase_rad":float(phase[yy,xx]),"dominant_kx":float(kx[yy,xx]),"dominant_ky":float(ky[yy,xx]),"dominant_power":float(power[yy,xx])})
    r=np.hypot(x-cx,y-cy); ap=r<=max(2.0,min(nx,ny)*0.04); ann=(r>=min(nx,ny)*0.10)&(r<=min(nx,ny)*0.22)
    bg=a[ann]; bgmed=float(np.median(bg)); mad=float(np.median(np.abs(bg-bgmed))); sig=1.4826*mad if mad>0 else float(np.std(bg))
    out.update({"center_mean":float(np.mean(a[ap])),"background_median":bgmed,"center_robust_z":float((np.mean(a[ap])-bgmed)/sig) if sig>0 else None})
    return out

def correlations(images):
    out={}; names=sorted(images)
    for i,l in enumerate(names):
        for r in names[i+1:]:
            a=clean(images[l]); b=clean(images[r])
            if a.shape!=b.shape: continue
            x=(a-a.mean()).ravel(); y=(b-b.mean()).ravel(); d=float(np.linalg.norm(x)*np.linalg.norm(y))
            out[f"{l}__{r}"]=float(np.dot(x,y)/d) if d else None
    return out

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("fits",nargs="+"); ap.add_argument("--output",required=True); ns=ap.parse_args()
    maps={}; images={}
    for raw in ns.fits:
        p=Path(raw); im,h=load_2d(p); images[p.stem]=im
        maps[p.stem]={"file":p.name,"bytes":p.stat().st_size,"sha256":sha256_file(p),
                      "header":{k:h.get(k) for k in ("OBJECT","TELESCOP","INSTRUME","BUNIT","CTYPE1","CTYPE2","CRVAL1","CRVAL2","CDELT1","CDELT2","NAXIS1","NAXIS2") if k in h},
                      "markers":markers(im)}
    rec={"schema":"PBH_CMB_PHASE_MARKER_RECEIPT_V0_1","status":"ANALYZED","target":{"ra_deg":107.576646,"dec_deg":-25.907579},
         "maps":maps,"pairwise_pixel_correlations":correlations(images),
         "interpretation":"DESCRIPTIVE_ONLY_NOT_PBH_DETECTION","authority":"CANDIDATE_ONLY","canon_allowed":False}
    Path(ns.output).parent.mkdir(parents=True,exist_ok=True); Path(ns.output).write_text(json.dumps(rec,indent=2,sort_keys=True)+"\n")
    print(json.dumps(rec,sort_keys=True))
if __name__=="__main__": main()
