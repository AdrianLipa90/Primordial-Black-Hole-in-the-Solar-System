#!/usr/bin/env python3
from __future__ import annotations
import json, math, random, statistics

OBLIQUITY_DEG=23.439291111
MSUN_TO_MEARTH=332946.0
A_CRIT_AU=170.0

# Stable/metastable TNOs with a>170 AU from Siraj, Chyba & Tremaine (2024),
# Table 1. Columns: name,a,q,i,varpi,Omega,survival,shift,metastable.
OBJECTS=[
("2016 TP120",170.7,40.4,32.64,116.9,126.69,82,11,False),
("2016 QV89",177.8,40.1,21.36,93.81,173.1,98,10,False),
("2015 UN105",194.2,41.4,37.02,1.11,129.32,98,9,False),
("2003 SS422",200.7,39.6,16.78,358.14,150.99,88,34,True),
("2013 UT15",208.7,44.1,10.62,83.5,191.91,100,0,False),
("2000 CR105",222.1,43.9,22.77,84.62,128.3,100,2,False),
("2012 VP113",273.9,80.6,24.02,24.92,90.85,100,0,False),
("2013 FT28",281.5,43.4,17.43,258.71,217.69,100,17,False),
("2018 VM35",295.3,44.7,8.47,135.6,192.23,87,13,False),
("2014 WB556",301.7,42.8,24.15,350.59,114.89,100,16,False),
("2014 SR349",324.7,47.5,17.94,16.23,35.02,100,1,False),
("2010 GB174",339.3,48.4,21.57,117.49,130.77,100,1,False),
("2004 VN112",355.3,47.3,25.5,33.05,66.04,100,4,False),
("2016 SD106",386.9,42.7,4.8,22.26,219.5,82,41,True),
("2015 RX245",467.2,45.7,12.1,73.26,8.58,100,33,True),
("2015 BP519",484.7,35.4,54.12,122.77,135.15,84,42,True),
("2013 RA109",527.1,46.0,12.4,7.91,104.82,100,39,True),
("2003 VB12",541.6,76.4,11.93,95.14,144.3,100,0,False),
("2021 DK18",735.2,44.5,15.45,196.98,322.22,82,26,True),
("2013 SY99",884.5,50.1,4.21,61.31,29.53,100,23,True),
]

GIANTS=[
("Jupiter",9.5458e-4,5.204),
("Saturn",2.858e-4,9.583),
("Uranus",4.366e-5,19.191),
("Neptune",5.151e-5,30.07),
]

def unit(v):
    n=math.sqrt(sum(x*x for x in v))
    return tuple(x/n for x in v)

def mean(vs):
    return tuple(sum(v[j] for v in vs)/len(vs) for j in range(3))

def ecc_vector(a,q,i_deg,varpi_deg,Omega_deg):
    e=1.0-q/a
    omega=math.radians((varpi_deg-Omega_deg)%360.0)
    O=math.radians(Omega_deg); i=math.radians(i_deg)
    x=e*(math.cos(O)*math.cos(omega)-math.sin(O)*math.sin(omega)*math.cos(i))
    y=e*(math.sin(O)*math.cos(omega)+math.cos(O)*math.sin(omega)*math.cos(i))
    z=e*(math.sin(omega)*math.sin(i))
    return (x,y,z),e

def ecl_to_eq(v):
    eps=math.radians(OBLIQUITY_DEG)
    x,y,z=unit(v)
    y2=y*math.cos(eps)-z*math.sin(eps)
    z2=y*math.sin(eps)+z*math.cos(eps)
    q=unit((x,y2,z2))
    ra=math.degrees(math.atan2(q[1],q[0]))%360.0
    dec=math.degrees(math.asin(q[2]))
    return ra,dec

def lonlat(v):
    x,y,z=unit(v)
    return math.degrees(math.atan2(y,x))%360.0,math.degrees(math.asin(z))

def angle(a,b):
    a=unit(a); b=unit(b)
    d=max(-1.0,min(1.0,sum(x*y for x,y in zip(a,b))))
    return math.degrees(math.acos(d))

def pct(vals,p):
    vals=sorted(vals)
    x=(len(vals)-1)*p
    lo=math.floor(x); hi=math.ceil(x)
    if lo==hi: return vals[lo]
    t=x-lo
    return vals[lo]*(1-t)+vals[hi]*t

def circular_mean(values,weights):
    x=sum(w*math.cos(math.radians(a)) for a,w in zip(values,weights))
    y=sum(w*math.sin(math.radians(a)) for a,w in zip(values,weights))
    return math.degrees(math.atan2(y,x))%360.0, math.hypot(x,y)/sum(weights)

def bootstrap(vs,nominal,n=20000):
    rng=random.Random(20260928)
    offsets=[]
    for _ in range(n):
        sample=[vs[rng.randrange(len(vs))] for _ in vs]
        m=mean(sample)
        if sum(x*x for x in m)<=0: continue
        anti=tuple(-x for x in unit(m))
        offsets.append(angle(anti,nominal))
    return {"n":len(offsets),"p50_deg":pct(offsets,0.50),"p68_deg":pct(offsets,0.68),"p95_deg":pct(offsets,0.95)}

def secular_strength_relation(a_crit=A_CRIT_AU):
    S=sum(m*a*a for _,m,a in GIANTS)
    mu_msun_au3=S/(a_crit**5)
    mu_mearth_au3=mu_msun_au3*MSUN_TO_MEARTH
    return S,mu_msun_au3,mu_mearth_au3

def mass_from_orbit(ap,ep,a_crit=A_CRIT_AU):
    S,_,_=secular_strength_relation(a_crit)
    mp_msun=S*(ap**3)*((1-ep*ep)**1.5)/(a_crit**5)
    return mp_msun*MSUN_TO_MEARTH

def main():
    ev=[]; varpi=[]; eccs=[]
    for name,a,q,i,vp,O,surv,shift,meta in OBJECTS:
        v,e=ecc_vector(a,q,i,vp,O)
        ev.append(v); varpi.append(vp); eccs.append(e)

    m=mean(ev)
    mean_e=math.sqrt(sum(x*x for x in m))
    anti=tuple(-x for x in unit(m))
    anti_lon,anti_lat=lonlat(anti)
    anti_ra,anti_dec=ecl_to_eq(anti)

    vp_e,R_e=circular_mean(varpi,eccs)
    vp_eq,R_eq=circular_mean(varpi,[1.0]*len(varpi))
    vp_surv,R_surv=circular_mean(varpi,[e*s/100.0 for e,(*_,s,sh,mt) in zip(eccs,OBJECTS)])

    stable=[row for row in OBJECTS if not row[-1]]
    st_vp=[r[4] for r in stable]
    st_e=[1-r[2]/r[1] for r in stable]
    vp_st,R_st=circular_mean(st_vp,st_e)

    S,mu_sun,mu_earth=secular_strength_relation()
    grid={}
    for ep in (0.0,0.2,0.29,0.4):
        grid[str(ep)]={str(ap):mass_from_orbit(ap,ep) for ap in (250,290,330,380)}

    out={
      "schema":"PBH_RELATIONAL_ORBITAL_INVERSION_V0_1",
      "status":"ANALYZED",
      "source_sample":{
        "n":len(OBJECTS),
        "criterion":"published stable/metastable TNOs with a>170 AU from Siraj-Chyba-Tremaine Table 1",
      },
      "reproduction":{
        "mean_eccentricity_vector":list(m),
        "mean_eccentricity_vector_magnitude":mean_e,
        "published_reference_value_approx":0.44,
      },
      "relational_axis":{
        "observed_e_weighted_varpi_deg":vp_e,
        "anti_aligned_varpi_deg":(vp_e+180.0)%360.0,
        "equal_weight_anti_varpi_deg":(vp_eq+180.0)%360.0,
        "survival_weighted_anti_varpi_deg":(vp_surv+180.0)%360.0,
        "stable_only_e_weighted_anti_varpi_deg":(vp_st+180.0)%360.0,
        "full_3d_anti_eccentricity_axis":{
          "ecliptic_lon_deg":anti_lon,
          "ecliptic_lat_deg":anti_lat,
          "ra_deg":anti_ra,
          "dec_deg":anti_dec,
        },
        "bootstrap_object_resampling":bootstrap(ev,anti),
        "interpretation":"anti-apsidal latent direction; not instantaneous PBH sky position"
      },
      "secular_scale_relation":{
        "a_crit_au":A_CRIT_AU,
        "giant_planet_quadrupole_sum_msun_au2":S,
        "mu_equiv_msun_per_au3":mu_sun,
        "mu_equiv_mearth_per_au3":mu_earth,
        "definition":"mu_equiv = m_p/[a_p^3 (1-e_p^2)^(3/2)] from leading-order precession-scale matching",
        "mass_grid_mearth":grid,
        "boundary":"leading-order secular scaling only; full N-body fit still required"
      },
      "canon_allowed":False
    }
    print(json.dumps(out,sort_keys=True,indent=2))

if __name__=="__main__":
    main()
