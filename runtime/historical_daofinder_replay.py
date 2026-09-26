#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import io
import json
import math
from urllib.parse import quote

import numpy as np
import requests
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.wcs import WCS
from photutils.detection import DAOStarFinder

SOURCE_REPO="AdrianLipa90/Metatime"
SOURCE_COMMIT="aa9fe46d2ca0a93e804e78ea843bbc296aeb591d"
RA_DEG=107.57664605896267
DEC_DEG=-25.907578906023748
TARGET=SkyCoord(ra=RA_DEG*u.deg,dec=DEC_DEG*u.deg,frame="icrs")
TARGET_GAL=TARGET.galactic

FILES={
 "planck_100_R3":("107_-26_Equ_J2000_planck_100Rel3.fits","1bef9697e314ad212d8aa2d1eec606968e189ead61a2d07e3d9eebec691409c5"),
 "planck_143_R1":("107_-26_Equ_J2000_planck_143Rel1.fits","893e21091e6c9b770c9acc17a88ae874bc57ccef8a0e3dd35c5bce5552de1f8b"),
 "planck_143_R2":("107_-26_Equ_J2000_planck_143Rel2.fits","63ebab928422961985d69b537744785d8321e5f023dfe523bb6586be524f1748"),
 "planck_143_R3":("107_-26_Equ_J2000_planck_143Rel3.fits","61770ad10ac6eb20b85a6587da21149b5e46d69b33273ff9f5eea93ed09d44a1"),
 "planck_217_R3":("107_-26_Equ_J2000_planck_217Rel3.fits","ba03fb5d2aaedfc56de6079529cd6bdc22ef875dd5019d6d2abcb7c5c0e7cd1a"),
}

def finite_image(data):
    a=np.asarray(data,dtype=float)
    med=float(np.nanmedian(a[np.isfinite(a)]))
    return np.where(np.isfinite(a),a,med)

def fetch(name):
    url=f"https://raw.githubusercontent.com/{SOURCE_REPO}/{SOURCE_COMMIT}/data/{quote(name)}"
    r=requests.get(url,timeout=120,headers={"User-Agent":"FERRET-PBH-historical-detector-replay/0.1"})
    if r.status_code!=200 or not r.content.startswith(b"SIMPLE"):
        raise RuntimeError(f"bad response HTTP {r.status_code}")
    return r.content,url

def serialize_sources(tab,limit=12):
    if tab is None:
        return []
    rows=[]
    for row in tab[:limit]:
        rows.append({
          "id":int(row["id"]),
          "xcentroid":float(row["xcentroid"]),
          "ycentroid":float(row["ycentroid"]),
          "peak":float(row["peak"]),
          "flux":float(row["flux"])
        })
    return rows

def main():
    receipt={
      "schema":"PBH_HISTORICAL_DAOSTARFINDER_REPLAY_V0_1",
      "status":"ANALYZED",
      "source_script":"AdrianLipa90/Metatime@aa9fe46d2ca0a93e804e78ea843bbc296aeb591d/simulations/PBH/primordialblackhole.py",
      "target":{"ra_deg":RA_DEG,"dec_deg":DEC_DEG,"gal_l_deg":float(TARGET_GAL.l.deg),"gal_b_deg":float(TARGET_GAL.b.deg)},
      "original_algorithm":{"sigma_clip":3.0,"dao_fwhm_pixels":3.0,"dao_threshold_sigma":5.0,"proximity_pixels":5.0},
      "maps":{},
      "interpretation":"ALGORITHM_REPLAY_ONLY_NOT_PBH_DETECTION",
      "authority":"CANDIDATE_ONLY",
      "canon_allowed":False,
    }
    for key,(name,expected_sha) in FILES.items():
        entry={"file":name,"expected_sha256":expected_sha}
        try:
            payload,url=fetch(name)
            sha=hashlib.sha256(payload).hexdigest()
            entry.update({"source_url":url,"bytes":len(payload),"sha256":sha,"sha256_match":sha==expected_sha})
            with fits.open(io.BytesIO(payload),memmap=False) as hdul:
                h=hdul[0].header
                data=finite_image(hdul[0].data)
            w=WCS(h)
            # Exact historical call: RA/Dec numeric values passed into GLON/GLAT WCS.
            wrong_x,wrong_y=w.all_world2pix(RA_DEG,DEC_DEG,1)
            # Correct physical target: transform ICRS target to Galactic first.
            right_x,right_y=w.all_world2pix(float(TARGET_GAL.l.deg),float(TARGET_GAL.b.deg),1)

            mean,median,std=sigma_clipped_stats(data,sigma=3.0)
            sub=data-median
            finder=DAOStarFinder(fwhm=3.0,threshold=5.0*std)
            sources=finder(sub)

            if sources is None or len(sources)==0:
                nearest_wrong=None
                nearest_right=None
                nsrc=0
                near_right_count=0
            else:
                dx=np.asarray(sources["xcentroid"],float)
                dy=np.asarray(sources["ycentroid"],float)
                wrong_dist=np.hypot(dx-float(wrong_x),dy-float(wrong_y))
                right_dist=np.hypot(dx-float(right_x),dy-float(right_y))
                nearest_wrong=float(np.nanmin(wrong_dist)) if np.isfinite(wrong_dist).any() else None
                nearest_right=float(np.min(right_dist))
                nsrc=int(len(sources))
                near_right_count=int(np.count_nonzero(right_dist<5.0))

            entry.update({
              "wcs":{"CTYPE1":h.get("CTYPE1"),"CTYPE2":h.get("CTYPE2"),"CRVAL1":h.get("CRVAL1"),"CRVAL2":h.get("CRVAL2"),"CDELT1":h.get("CDELT1"),"CDELT2":h.get("CDELT2"),"BUNIT":h.get("BUNIT")},
              "historical_untransformed_target_pixel_origin1":[
                  float(wrong_x) if np.isfinite(wrong_x) else None,
                  float(wrong_y) if np.isfinite(wrong_y) else None
              ],
              "historical_target_pixel_is_finite":bool(np.isfinite(wrong_x) and np.isfinite(wrong_y)),
              "correct_galactic_target_pixel_origin1":[float(right_x),float(right_y)],
              "sigma_clipped":{"mean":float(mean),"median":float(median),"std":float(std),"threshold_5sigma":float(5.0*std)},
              "dao_source_count":nsrc,
              "nearest_source_distance_pixels_wrong_target":nearest_wrong,
              "nearest_source_distance_pixels_correct_target":nearest_right,
              "sources_within_5px_correct_target":near_right_count,
              "historical_alert_if_literal_code":bool(nearest_wrong is not None and nearest_wrong<5.0),
              "corrected_alert_if_same_5px_rule":bool(nearest_right is not None and nearest_right<5.0),
              "brightest_sources_preview":serialize_sources(sources),
            })
        except Exception as exc:
            entry["error"]=f"{type(exc).__name__}: {exc}"
        receipt["maps"][key]=entry

    receipt["bestiary"]={
      "OWL":"The original source detector passes ICRS RA/Dec numbers directly to a Galactic GLON/GLAT WCS, so its PBH target pixel is invalid or wrong.",
      "HOUND":"Replay separately measures the literal-code target and the correctly transformed Galactic target.",
      "BAT":"DAOStarFinder detections are generic local extrema in diffuse CMB frequency maps, not a PBH-specific observable.",
      "BELZEBUB":"A 5-pixel proximity alert cannot be promoted to PBH evidence without a calibrated null model and an observable tied to lensing/astrometry.",
      "GREMLIN":"HISTORICAL_POINT_SOURCE_TEST_REPLAYED_WITH_WCS_FIREWALL"
    }
    print("PBH_DAO_REPLAY="+json.dumps(receipt,sort_keys=True,separators=(",",":"),allow_nan=False),flush=True)

if __name__=="__main__":
    main()
