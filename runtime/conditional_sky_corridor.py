#!/usr/bin/env python3
from __future__ import annotations
import json, math
from datetime import datetime, timezone
import numpy as np
from astropy.time import Time
from astropy.coordinates import get_body_barycentric_posvel, solar_system_ephemeris
import astropy.units as u

AU_M=149597870700.0
DAY_S=86400.0
YEAR_D=365.25
GM_SUN=1.32712440018e20
MU=GM_SUN*DAY_S**2/AU_M**3  # AU^3/day^2
EPOCH=Time("2025-03-07T00:00:00",scale="utc")
DATES=[
    "2025-03-07T00:00:00",
    "2026-09-26T00:00:00",
    "2027-03-07T00:00:00",
    "2027-09-26T00:00:00",
]
R0=np.array([-158.399,500.044,-254.785],dtype=float)
VEL_AUYR={
 "A_OUTBOUND":np.array([-0.1766790385976647,0.09911727381056191,-0.050502613021396095]),
 "A_INBOUND":np.array([-0.09244534800482046,-0.1667969780395646,0.08498738917367973]),
 "B_OUTBOUND":np.array([-0.1716206865269884,0.08314874747624625,-0.08498738917367979]),
 "B_INBOUND":np.array([-0.08738699593414413,-0.18276550437388026,0.05050261302139604]),
}

def accel(r):
    d=float(np.linalg.norm(r))
    return -MU*r/(d**3)

def propagate(r0,v0_auyr,delta_days,step_days=0.25):
    r=r0.copy()
    v=v0_auyr/YEAR_D
    sign=1.0 if delta_days>=0 else -1.0
    rem=abs(float(delta_days))
    while rem>1e-12:
        h=sign*min(step_days,rem)
        a0=accel(r)
        r_new=r+v*h+0.5*a0*h*h
        a1=accel(r_new)
        v_new=v+0.5*(a0+a1)*h
        r,v=r_new,v_new
        rem-=abs(h)
    return r,v*YEAR_D

def radec(vec):
    x,y,z=map(float,vec)
    rr=math.sqrt(x*x+y*y+z*z)
    ra=math.degrees(math.atan2(y,x))%360.0
    dec=math.degrees(math.asin(z/rr))
    return ra,dec,rr

def sep_arcmin(ra1,dec1,ra2,dec2):
    a1,d1,a2,d2=map(math.radians,[ra1,dec1,ra2,dec2])
    c=math.sin(d1)*math.sin(d2)+math.cos(d1)*math.cos(d2)*math.cos(a1-a2)
    c=max(-1.0,min(1.0,c))
    return math.degrees(math.acos(c))*60.0

def earth_heliocentric_icrs(t):
    with solar_system_ephemeris.set("builtin"):
        ep,_=get_body_barycentric_posvel("earth",t)
        sp,_=get_body_barycentric_posvel("sun",t)
    v=(ep.xyz-sp.xyz).to_value(u.AU)
    return np.array(v,dtype=float)

def main():
    rows={}
    for name,v0 in VEL_AUYR.items():
        rows[name]={}
        for ds in DATES:
            t=Time(ds,scale="utc")
            dt=(t-EPOCH).to_value(u.day)
            r,v=propagate(R0,v0,dt)
            era=earth_heliocentric_icrs(t)
            g=r-era
            hra,hdec,hdist=radec(r)
            gra,gdec,gdist=radec(g)
            rows[name][ds]={
              "heliocentric_position_AU":r.tolist(),
              "heliocentric_velocity_AU_yr":v.tolist(),
              "heliocentric_ra_deg":hra,
              "heliocentric_dec_deg":hdec,
              "heliocentric_distance_AU":hdist,
              "geocentric_ra_deg":gra,
              "geocentric_dec_deg":gdec,
              "geocentric_distance_AU":gdist,
              "earth_heliocentric_icrs_AU":era.tolist(),
              "parallax_offset_arcmin":sep_arcmin(hra,hdec,gra,gdec),
            }

    envelope={}
    for ds in DATES:
        vals=[(n,rows[n][ds]["geocentric_ra_deg"],rows[n][ds]["geocentric_dec_deg"]) for n in rows]
        maxsep=0.0; pair=None
        for i in range(len(vals)):
            for j in range(i+1,len(vals)):
                s=sep_arcmin(vals[i][1],vals[i][2],vals[j][1],vals[j][2])
                if s>maxsep:
                    maxsep=s; pair=[vals[i][0],vals[j][0]]
        ras=np.array([v[1] for v in vals]); decs=np.array([v[2] for v in vals])
        envelope[ds]={
          "max_branch_pair_separation_arcmin":maxsep,
          "max_separation_pair":pair,
          "ra_deg_min":float(ras.min()),"ra_deg_max":float(ras.max()),
          "dec_deg_min":float(decs.min()),"dec_deg_max":float(decs.max()),
        }

    receipt={
      "schema":"PBH_CONDITIONAL_SKY_CORRIDOR_V0_1",
      "status":"ANALYZED",
      "assumption":{
        "published_cartesian_vector_epoch":"2025-03-07T00:00:00 UTC",
        "epoch_status":"PROVISIONAL_NOT_RECOVERED_FROM_ORIGINAL_SOURCE",
        "frame_assumption":"published Cartesian vector interpreted in ICRS-equatorial axes because its direct atan2/asin conversion reproduces the published RA/Dec",
        "dynamics":"Sun-only two-body propagation; Earth position from Astropy builtin ephemeris for observer parallax",
        "branch_policy":"all four exact Kepler-consistent velocity/orientation branches retained"
      },
      "branches":rows,
      "envelope":envelope,
      "interpretation":"CONDITIONAL_SEARCH_CORRIDOR_ONLY",
      "hard_stop":"Do not call any row a current PBH position unless the original simulation epoch/state provenance is recovered or independently justified.",
      "authority":"CANDIDATE_ONLY",
      "canon_allowed":False
    }
    print("PBH_SKY_CORRIDOR="+json.dumps(receipt,sort_keys=True,separators=(",",":"),allow_nan=False),flush=True)

if __name__=="__main__":
    main()
