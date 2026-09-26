#!/usr/bin/env python3
from __future__ import annotations
import cmath, hashlib, io, json, math, sys, time
from urllib.parse import urlencode
import numpy as np
import requests
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u

TARGET_L=238.1709218507047
TARGET_B=-7.665782897020605
TAU=2*math.pi
WIDTH=128
FOV_DEG=8.0
PIX_ARCMIN=FOV_DEG*60/WIDTH
COMMON_FWHM=40.0
ENDPOINT="https://alasky.cds.unistra.fr/hips-image-services/hips2fits"
SOURCES={
 "planck_smica":("CDS/P/PLANCK/R3/CMB",5.0),
 "planck_100":("CDS/P/PLANCK/R3/HFI100",9.66),
 "planck_143":("CDS/P/PLANCK/R3/HFI143",7.27),
 "planck_217":("CDS/P/PLANCK/R3/HFI217",5.01),
 "wmap_V":("CDS/P/WMAP/V/9yr",21.0),
 "wmap_W":("CDS/P/WMAP/W/9yr",13.2),
}
FEATURES=("phase_coherence","phase_entropy","A1","A2","A3","A6")


def finite(a):
    a=np.asarray(a,float); m=np.isfinite(a)
    if not m.any(): raise RuntimeError("no finite pixels")
    med=float(np.nanmedian(a[m]))
    return np.where(m,a,med)


def smooth(a,native):
    a=finite(a)
    if native>=COMMON_FWHM:return a
    extra=math.sqrt(COMMON_FWHM**2-native**2)
    sigma=extra/(2.354820045*PIX_ARCMIN)
    ny,nx=a.shape
    fy=np.fft.fftfreq(ny); fx=np.fft.fftfreq(nx)
    ky,kx=np.meshgrid(fy,fx,indexing="ij")
    tr=np.exp(-2*(math.pi**2)*(sigma**2)*(kx*kx+ky*ky))
    return np.fft.ifft2(np.fft.fft2(a)*tr).real


def marker(a):
    a=finite(a); ny,nx=a.shape
    c=a-float(np.median(a)); F=np.fft.fftshift(np.fft.fft2(c*np.outer(np.hanning(ny),np.hanning(nx))))
    amp=np.abs(F); pwr=amp*amp; ph=np.angle(F)
    yy,xx=np.indices(a.shape); cy,cx=(ny-1)/2,(nx-1)/2
    kx=(xx-cx)/nx; ky=(yy-cy)/ny; kr=np.hypot(kx,ky); th=np.arctan2(ky,kx)
    band=(kr>0.03)&(kr<0.20); aw=amp[band]; pw=pwr[band]
    coh=float(abs(np.sum(aw*np.exp(1j*ph[band])))/np.sum(aw))
    hist,_=np.histogram(ph[band],bins=24,range=(-np.pi,np.pi),weights=aw)
    q=hist/hist.sum(); q=q[q>0]
    out={"phase_coherence":coh,"phase_entropy":float(-np.sum(q*np.log(q))/np.log(24))}
    for n in (1,2,3,6):
        out[f"A{n}"]=float(abs(np.sum(pw*np.exp(1j*n*th[band])))/np.sum(pw))
    return out


def center_z(a):
    a=finite(a); ny,nx=a.shape; yy,xx=np.indices(a.shape); cy,cx=(ny-1)/2,(nx-1)/2
    r=np.hypot(xx-cx,yy-cy)*(FOV_DEG/WIDTH)
    ap=r<=0.25; ann=(r>=0.75)&(r<=1.50)
    bg=a[ann]; med=float(np.median(bg)); mad=float(np.median(np.abs(bg-med)))
    sig=1.4826*mad if mad>0 else float(np.std(bg))
    return float((np.mean(a[ap])-med)/sig) if sig else None


def corr(a,b):
    x=(finite(a)-np.mean(a)).ravel(); y=(finite(b)-np.mean(b)).ravel()
    d=float(np.linalg.norm(x)*np.linalg.norm(y))
    return float(np.dot(x,y)/d) if d else None


def t36(m):
    base=[TAU*(float(m[k])%1.0) for k in FEATURES]
    return tuple(base*6)


def delta(a,b):
    return (a-b+math.pi)%TAU-math.pi


def phase_coh(a,b):
    z=sum(cmath.exp(1j*delta(x,y)) for x,y in zip(a,b))/36
    return abs(z)


def fetch(hips,ra,dec):
    params={"hips":hips,"width":WIDTH,"height":WIDTH,"projection":"TAN","fov":FOV_DEG,
            "ra":ra,"dec":dec,"coordsys":"icrs","rotation_angle":0,"format":"fits"}
    url=ENDPOINT+"?"+urlencode(params)
    r=requests.get(url,timeout=120,headers={"User-Agent":"GREMLIN-FERRET-PBH-NULL/0.1 read-only"})
    if r.status_code!=200 or not r.content.startswith(b"SIMPLE"):
        raise RuntimeError(f"bad FITS response {r.status_code}")
    with fits.open(io.BytesIO(r.content),memmap=False) as hdul:
        a=np.asarray(hdul[0].data,float)
        while a.ndim>2:a=a[0]
    return finite(a),hashlib.sha256(r.content).hexdigest()


def analyze_position(l_deg,b_deg):
    sky=SkyCoord(l=l_deg*u.deg,b=b_deg*u.deg,frame="galactic").icrs
    ra=float(sky.ra.deg); dec=float(sky.dec.deg)
    imgs={}; marks={}; z={}; hashes={}
    for name,(hips,fwhm) in SOURCES.items():
        a,h=fetch(hips,ra,dec); s=smooth(a,fwhm)
        q=WIDTH//4; central=s[q:3*q,q:3*q]
        imgs[name]=s; marks[name]=marker(central); z[name]=center_z(s); hashes[name]=h
    t={n:t36(m) for n,m in marks.items()}
    core=["planck_smica","planck_100","planck_143","wmap_V","wmap_W"]
    pcs=[]
    for i,a in enumerate(core):
        for b in core[i+1:]: pcs.append(phase_coh(t[a],t[b]))
    return {
      "l_deg":l_deg,"b_deg":b_deg,"ra_deg":ra,"dec_deg":dec,
      "center_z":z,
      "marker":marks,
      "gremlin_core_mean_pairwise_coherence":float(np.mean(pcs)),
      "smica_wmapV_image_corr":corr(imgs["planck_smica"],imgs["wmap_V"]),
      "smica_wmapW_image_corr":corr(imgs["planck_smica"],imgs["wmap_W"]),
      "wmapVW_image_corr":corr(imgs["wmap_V"],imgs["wmap_W"]),
      "hashes":hashes,
    }


def empirical_upper(target,controls):
    return (1+sum(x>=target for x in controls))/(1+len(controls))


def empirical_abs(target,controls):
    return (1+sum(abs(x)>=abs(target) for x in controls))/(1+len(controls))


def main():
    longitudes=[float(x) for x in range(0,360,15) if abs(((x-TARGET_L+180)%360)-180)>=30]
    target=analyze_position(TARGET_L,TARGET_B)
    controls=[]
    for i,l in enumerate(longitudes):
        controls.append(analyze_position(l,TARGET_B))
        print(f"CONTROL {i+1}/{len(longitudes)} l={l}",flush=True)
        time.sleep(0.05)
    out={
      "schema":"PBH_CMB_MATCHED_LATITUDE_NULL_V0_1",
      "status":"ANALYZED",
      "target":target,
      "control_definition":{"gal_b_deg":TARGET_B,"longitude_grid_deg":longitudes,"exclusion_from_target_deg":30.0},
      "n_controls":len(controls),
      "empirical":{
        "abs_smica_center_z":empirical_abs(target["center_z"]["planck_smica"],[c["center_z"]["planck_smica"] for c in controls]),
        "gremlin_core_mean_pairwise_coherence_upper":empirical_upper(target["gremlin_core_mean_pairwise_coherence"],[c["gremlin_core_mean_pairwise_coherence"] for c in controls]),
        "smica_wmapV_image_corr_upper":empirical_upper(target["smica_wmapV_image_corr"],[c["smica_wmapV_image_corr"] for c in controls]),
        "smica_wmapW_image_corr_upper":empirical_upper(target["smica_wmapW_image_corr"],[c["smica_wmapW_image_corr"] for c in controls]),
      },
      "controls_summary":[
        {"l_deg":c["l_deg"],"ra_deg":c["ra_deg"],"dec_deg":c["dec_deg"],
         "smica_center_z":c["center_z"]["planck_smica"],
         "gremlin_core_mean_pairwise_coherence":c["gremlin_core_mean_pairwise_coherence"],
         "smica_wmapV_image_corr":c["smica_wmapV_image_corr"],
         "smica_wmapW_image_corr":c["smica_wmapW_image_corr"]}
        for c in controls
      ],
      "interpretation":"EMPIRICAL_MATCHED_LATITUDE_CONTROL_ONLY_NOT_FORMAL_SURVEY_PVALUE",
      "canon_allowed":False,
    }
    with open("/tmp/null_receipt.json","w") as f: json.dump(out,f,sort_keys=True,indent=2)
    print("PBH_NULL_RECEIPT="+json.dumps(out,sort_keys=True,separators=(",",":")),flush=True)
    return 0

if __name__=="__main__":sys.exit(main())
