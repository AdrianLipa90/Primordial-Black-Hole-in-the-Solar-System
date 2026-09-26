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
WIDTH=256
FOV_DEG=8.0
PIX_ARCMIN=FOV_DEG*60.0/WIDTH
COMMON_FWHM=40.0
ENDPOINT="https://alasky.cds.unistra.fr/hips-image-services/hips2fits"
SOURCES={
 "planck_smica":("CDS/P/PLANCK/R3/CMB","CMB_COMPONENT",5.0),
 "planck_100":("CDS/P/PLANCK/R3/HFI100","CMB_FREQUENCY",9.66),
 "planck_143":("CDS/P/PLANCK/R3/HFI143","CMB_FREQUENCY",7.27),
 "planck_217":("CDS/P/PLANCK/R3/HFI217","CMB_FREQUENCY",5.01),
 "wmap_V":("CDS/P/WMAP/V/9yr","CMB_FREQUENCY",21.0),
 "wmap_W":("CDS/P/WMAP/W/9yr","CMB_FREQUENCY",13.2),
}

def finite(a):
    a=np.asarray(a,float); m=np.isfinite(a)
    if not m.any(): raise RuntimeError("no finite pixels")
    return np.where(m,a,float(np.nanmedian(a[m])))

def smooth(a,native):
    a=finite(a)
    if native>=COMMON_FWHM:return a
    extra=math.sqrt(COMMON_FWHM**2-native**2)
    sigma=extra/(2.354820045*PIX_ARCMIN)
    ny,nx=a.shape; fy=np.fft.fftfreq(ny); fx=np.fft.fftfreq(nx)
    ky,kx=np.meshgrid(fy,fx,indexing="ij")
    transfer=np.exp(-2.0*(math.pi**2)*(sigma**2)*(kx*kx+ky*ky))
    return np.fft.ifft2(np.fft.fft2(a)*transfer).real

def phase36_sector_cell(image):
    a=finite(image); ny,nx=a.shape
    c=a-float(np.median(a)); F=np.fft.fftshift(np.fft.fft2(c*np.outer(np.hanning(ny),np.hanning(nx))))
    amp=np.abs(F); ph=np.angle(F)
    yy,xx=np.indices(a.shape); cy,cx=(ny-1)/2,(nx-1)/2
    kx=(xx-cx)/nx; ky=(yy-cy)/ny; kr=np.hypot(kx,ky); th=np.arctan2(ky,kx)
    re=np.linspace(0.02,0.20,7); ae=np.linspace(-math.pi,math.pi,7)
    phases=[]; coherence=[]
    for ri in range(6):
        for ai in range(6):
            mask=(kr>=re[ri])&(kr<re[ri+1])&(th>=ae[ai])&(th<ae[ai+1])
            w=amp[mask]; p=ph[mask]
            if w.size==0 or float(w.sum())<=1e-18:
                phases.append(0.0); coherence.append(0.0); continue
            z=np.sum(w*np.exp(1j*p))
            phases.append(float(np.angle(z)%(2*math.pi)))
            coherence.append(float(abs(z)/w.sum()))
    return {"phase36":phases,"global_sector_coherence":float(np.mean(coherence))}

def circular_delta(a,b):
    return math.atan2(math.sin(a-b),math.cos(a-b))

def moire_overlap(a,b):
    z=sum(cmath.exp(1j*circular_delta(x,y)) for x,y in zip(a,b))/36.0
    return abs(z)

def phase_rms(a,b):
    return math.sqrt(sum(circular_delta(x,y)**2 for x,y in zip(a,b))/36.0)

def project_hologram(cells):
    amps=[]
    for j in range(36):
        z=sum(cmath.exp(1j*c[j]) for c in cells)
        amps.append(abs(z)/len(cells))
    return float(np.mean(amps))

def center_z(a):
    a=finite(a); ny,nx=a.shape; yy,xx=np.indices(a.shape); cy,cx=(ny-1)/2,(nx-1)/2
    r=np.hypot(xx-cx,yy-cy)*(FOV_DEG/WIDTH)
    ap=r<=0.25; ann=(r>=0.75)&(r<=1.50)
    bg=a[ann]; med=float(np.median(bg)); mad=float(np.median(np.abs(bg-med)))
    sig=1.4826*mad if mad>0 else float(np.std(bg))
    return float((np.mean(a[ap])-med)/sig) if sig else None

def fetch(hips,ra,dec):
    params={"hips":hips,"width":WIDTH,"height":WIDTH,"projection":"TAN","fov":FOV_DEG,
            "ra":ra,"dec":dec,"coordsys":"icrs","rotation_angle":0,"format":"fits"}
    url=ENDPOINT+"?"+urlencode(params)
    r=requests.get(url,timeout=120,headers={"User-Agent":"GREMLIN-FERRET-PBH-EXACT-NULL/0.1 read-only"})
    if r.status_code!=200 or not r.content.startswith(b"SIMPLE"):
        raise RuntimeError(f"bad FITS response {r.status_code}")
    with fits.open(io.BytesIO(r.content),memmap=False) as hdul:
        a=np.asarray(hdul[0].data,float)
        while a.ndim>2:a=a[0]
    return finite(a),hashlib.sha256(r.content).hexdigest()

def analyze(l_deg,b_deg):
    sky=SkyCoord(l=l_deg*u.deg,b=b_deg*u.deg,frame="galactic").icrs
    ra,dec=float(sky.ra.deg),float(sky.dec.deg)
    cells={}; center={}; hashes={}
    for name,(hips,kind,fwhm) in SOURCES.items():
        a,h=fetch(hips,ra,dec); s=smooth(a,fwhm); q=WIDTH//4; central=s[q:3*q,q:3*q]
        cells[name]=phase36_sector_cell(central)
        center[name]=center_z(s); hashes[name]=h
    freq=[name for name,(_,kind,_) in SOURCES.items() if kind=="CMB_FREQUENCY"]
    all_names=sorted(cells)
    pair={}
    for i,a in enumerate(all_names):
        for b in all_names[i+1:]:
            pair[f"{a}__{b}"]={"moire_overlap":moire_overlap(cells[a]["phase36"],cells[b]["phase36"]),
                                "phase_rms_distance_rad":phase_rms(cells[a]["phase36"],cells[b]["phase36"])}
    return {
      "l_deg":l_deg,"b_deg":b_deg,"ra_deg":ra,"dec_deg":dec,
      "center_z":center,
      "global_sector_coherence":{n:cells[n]["global_sector_coherence"] for n in cells},
      "frequency_hologram_global_coherence":project_hologram([cells[n]["phase36"] for n in freq]),
      "all_hologram_global_coherence":project_hologram([cells[n]["phase36"] for n in all_names]),
      "smica_wmapV_moire":pair["planck_smica__wmap_V"]["moire_overlap"],
      "smica_wmapW_moire":pair["planck_smica__wmap_W"]["moire_overlap"],
      "wmapVW_moire":pair["wmap_V__wmap_W"]["moire_overlap"],
      "pairwise":pair,"hashes":hashes,
    }

def tail_upper(target,controls):
    return (1+sum(x>=target for x in controls))/(1+len(controls))

def tail_abs(target,controls):
    return (1+sum(abs(x)>=abs(target) for x in controls))/(1+len(controls))

def main():
    ls=[float(x) for x in range(0,360,45) if abs(((x-TARGET_L+180)%360)-180)>=30]
    target=analyze(TARGET_L,TARGET_B)
    controls=[]
    for i,l in enumerate(ls):
        controls.append(analyze(l,TARGET_B))
        print(f"EXACT_CONTROL {i+1}/{len(ls)} l={l}",flush=True)
        time.sleep(0.05)
    out={
      "schema":"PBH_CMB_GREMLIN_T36_MATCHED_LATITUDE_PILOT_NULL_V0_1",
      "status":"ANALYZED",
      "gremlin_algebra":{"phase_cell":"RADIAL6_X_ANGULAR6_ROW_MAJOR","distance":"phase_rms_distance","overlap":"moire_overlap","projection":"project_hologram","source":"AdrianLipa90/GREMLIN gremlin_mcp/holographic_orbital_memory.py"},
      "target":target,
      "controls":[{"l_deg":c["l_deg"],"ra_deg":c["ra_deg"],"dec_deg":c["dec_deg"],
        "smica_center_z":c["center_z"]["planck_smica"],
        "frequency_hologram_global_coherence":c["frequency_hologram_global_coherence"],
        "all_hologram_global_coherence":c["all_hologram_global_coherence"],
        "smica_wmapV_moire":c["smica_wmapV_moire"],
        "smica_wmapW_moire":c["smica_wmapW_moire"],
        "wmapVW_moire":c["wmapVW_moire"]} for c in controls],
      "empirical":{
        "abs_smica_center_z":tail_abs(target["center_z"]["planck_smica"],[c["center_z"]["planck_smica"] for c in controls]),
        "frequency_hologram_coherence_upper":tail_upper(target["frequency_hologram_global_coherence"],[c["frequency_hologram_global_coherence"] for c in controls]),
        "all_hologram_coherence_upper":tail_upper(target["all_hologram_global_coherence"],[c["all_hologram_global_coherence"] for c in controls]),
        "smica_wmapV_moire_upper":tail_upper(target["smica_wmapV_moire"],[c["smica_wmapV_moire"] for c in controls]),
        "smica_wmapW_moire_upper":tail_upper(target["smica_wmapW_moire"],[c["smica_wmapW_moire"] for c in controls]),
      },
      "interpretation":"PILOT_MATCHED_LATITUDE_EXACT_GREMLIN_T36_NOT_FORMAL_SURVEY_PVALUE",
      "authority":"CANDIDATE_ONLY","canon_allowed":False,
    }
    with open("/tmp/exact_null_receipt.json","w") as f:json.dump(out,f,sort_keys=True,indent=2)
    print("PBH_EXACT_NULL_RECEIPT="+json.dumps(out,sort_keys=True,separators=(",",":")),flush=True)
    return 0

if __name__=="__main__":sys.exit(main())
