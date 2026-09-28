#!/usr/bin/env python3
from __future__ import annotations

import json
import math
import random
import statistics
import urllib.parse
import urllib.request
from typing import Any

JPL_QUERY = "https://ssd-api.jpl.nasa.gov/sbdb_query.api"
OBLIQUITY_DEG = 23.439291111
SEED = 20260928

FIELDS = [
    "full_name","epoch","a","e","q","i","om","w","ma",
    "sigma_a","sigma_e","sigma_q","sigma_i","sigma_om","sigma_w",
    "condition_code","data_arc","n_obs_used"
]

SAMPLE_DEFS = {
    "A150_Q30": {"a_min":150.0, "q_min":30.0},
    "A250_Q30": {"a_min":250.0, "q_min":30.0},
    "A250_Q35": {"a_min":250.0, "q_min":35.0},
    "A300_Q35": {"a_min":300.0, "q_min":35.0},
}

def fnum(v: Any) -> float | None:
    if v in (None, "", "null"):
        return None
    try:
        x=float(v)
    except Exception:
        return None
    return x if math.isfinite(x) else None

def inum(v: Any) -> int | None:
    x=fnum(v)
    return None if x is None else int(x)

def deg2rad(x: float) -> float:
    return math.radians(x)

def unit(v):
    n=math.sqrt(sum(x*x for x in v))
    if n <= 0:
        return None
    return tuple(x/n for x in v)

def mean_unit(vectors):
    s=[sum(v[j] for v in vectors) for j in range(3)]
    m=unit(s)
    if m is None:
        return None,0.0
    r=math.sqrt(sum(x*x for x in s))/len(vectors)
    return m,r

def dot(a,b):
    return sum(x*y for x,y in zip(a,b))

def angle_deg(a,b):
    d=max(-1.0,min(1.0,dot(a,b)))
    return math.degrees(math.acos(d))

def perihelion_vector(i_deg, om_deg, w_deg):
    i=deg2rad(i_deg); Om=deg2rad(om_deg); w=deg2rad(w_deg)
    # Standard heliocentric ecliptic direction of radius vector at true anomaly f=0.
    x=math.cos(Om)*math.cos(w)-math.sin(Om)*math.sin(w)*math.cos(i)
    y=math.sin(Om)*math.cos(w)+math.cos(Om)*math.sin(w)*math.cos(i)
    z=math.sin(w)*math.sin(i)
    return unit((x,y,z))

def pole_vector(i_deg, om_deg):
    i=deg2rad(i_deg); Om=deg2rad(om_deg)
    # Orbit angular-momentum unit vector in J2000 ecliptic coordinates.
    return unit((math.sin(i)*math.sin(Om), -math.sin(i)*math.cos(Om), math.cos(i)))

def ecl_lonlat(v):
    x,y,z=v
    lon=math.degrees(math.atan2(y,x))%360.0
    lat=math.degrees(math.asin(max(-1.0,min(1.0,z))))
    return lon,lat

def ecl_to_equatorial(v):
    eps=deg2rad(OBLIQUITY_DEG)
    x,y,z=v
    xe=x
    ye=y*math.cos(eps)-z*math.sin(eps)
    ze=y*math.sin(eps)+z*math.cos(eps)
    q=unit((xe,ye,ze))
    x,y,z=q
    ra=math.degrees(math.atan2(y,x))%360.0
    dec=math.degrees(math.asin(max(-1.0,min(1.0,z))))
    return ra,dec

def antipode(v):
    return tuple(-x for x in v)

def circular_mean_deg(values):
    xs=sum(math.cos(deg2rad(v)) for v in values)
    ys=sum(math.sin(deg2rad(v)) for v in values)
    ang=math.degrees(math.atan2(ys,xs))%360.0
    r=math.hypot(xs,ys)/len(values)
    return ang,r

def percentile(vals,p):
    if not vals:
        return None
    a=sorted(vals)
    if len(a)==1:
        return a[0]
    pos=(len(a)-1)*p
    lo=int(math.floor(pos)); hi=int(math.ceil(pos))
    if lo==hi:
        return a[lo]
    t=pos-lo
    return a[lo]*(1-t)+a[hi]*t

def bootstrap_axis(vectors, nominal_antipode, n=4000):
    rng=random.Random(SEED)
    offsets=[]
    ras=[]
    decs=[]
    for _ in range(n):
        sample=[vectors[rng.randrange(len(vectors))] for _ in vectors]
        m,r=mean_unit(sample)
        if m is None:
            continue
        a=antipode(m)
        offsets.append(angle_deg(a,nominal_antipode))
        ra,dec=ecl_to_equatorial(a)
        ras.append(ra); decs.append(dec)
    return {
        "n":len(offsets),
        "angular_offset_deg_p50":percentile(offsets,0.50),
        "angular_offset_deg_p68":percentile(offsets,0.68),
        "angular_offset_deg_p95":percentile(offsets,0.95),
        "note":"resampling spread only; does not include survey-selection bias or full orbital covariance"
    }

def leave_one_out(vectors, names, nominal_antipode):
    rows=[]
    if len(vectors)<3:
        return rows
    for k,name in enumerate(names):
        sub=[v for j,v in enumerate(vectors) if j!=k]
        m,r=mean_unit(sub)
        if m is None:
            continue
        a=antipode(m)
        ra,dec=ecl_to_equatorial(a)
        rows.append({
            "removed":name,
            "axis_shift_deg":angle_deg(a,nominal_antipode),
            "anti_axis_ra_deg":ra,
            "anti_axis_dec_deg":dec,
            "resultant":r
        })
    return rows

def fetch_rows():
    cdata=json.dumps({"AND":["a|GT|150","q|GT|30"]},separators=(",",":"))
    params={
        "fields":",".join(FIELDS),
        "sb-class":"TNO",
        "sb-cdata":cdata,
        "limit":"500"
    }
    url=JPL_QUERY+"?"+urllib.parse.urlencode(params)
    req=urllib.request.Request(url,headers={"User-Agent":"PBH-relational-orbit-inference/0.1"})
    with urllib.request.urlopen(req,timeout=60) as resp:
        payload=json.loads(resp.read().decode("utf-8"))
    fields=payload["fields"]
    out=[]
    for row in payload["data"]:
        d=dict(zip(fields,row))
        a=fnum(d.get("a")); q=fnum(d.get("q")); e=fnum(d.get("e"))
        inc=fnum(d.get("i")); om=fnum(d.get("om")); w=fnum(d.get("w"))
        if None in (a,q,e,inc,om,w):
            continue
        d2={
            "name":str(d.get("full_name")).strip(),
            "epoch":fnum(d.get("epoch")),
            "a":a,"e":e,"q":q,"i":inc,"om":om,"w":w,
            "ma":fnum(d.get("ma")),
            "sigma_a":fnum(d.get("sigma_a")),
            "sigma_e":fnum(d.get("sigma_e")),
            "sigma_q":fnum(d.get("sigma_q")),
            "sigma_i":fnum(d.get("sigma_i")),
            "sigma_om":fnum(d.get("sigma_om")),
            "sigma_w":fnum(d.get("sigma_w")),
            "condition_code":inum(d.get("condition_code")),
            "data_arc_days":fnum(d.get("data_arc")),
            "n_obs":inum(d.get("n_obs_used")),
        }
        d2["varpi_deg"]=(om+w)%360.0
        d2["perihelion_vector_ecl"]=perihelion_vector(inc,om,w)
        d2["pole_vector_ecl"]=pole_vector(inc,om)
        out.append(d2)
    return payload.get("signature"),url,out

def quality_ok(o):
    # Conservative reproducibility gate, not a claim that rejected objects are physically invalid.
    if o["data_arc_days"] is not None and o["data_arc_days"] < 365.0:
        return False
    if o["n_obs"] is not None and o["n_obs"] < 10:
        return False
    if o["condition_code"] is not None and o["condition_code"] > 6:
        return False
    return True

def analyze_sample(rows, spec):
    sub=[o for o in rows if o["a"]>spec["a_min"] and o["q"]>spec["q_min"] and quality_ok(o)]
    if len(sub)<3:
        return {"status":"INSUFFICIENT","n":len(sub),"objects":[o["name"] for o in sub]}

    pvec=[o["perihelion_vector_ecl"] for o in sub]
    poles=[o["pole_vector_ecl"] for o in sub]
    names=[o["name"] for o in sub]

    pmean,pR=mean_unit(pvec)
    polemean,poleR=mean_unit(poles)
    anti=antipode(pmean)

    peri_lon,peri_lat=ecl_lonlat(pmean)
    anti_lon,anti_lat=ecl_lonlat(anti)
    anti_ra,anti_dec=ecl_to_equatorial(anti)
    pole_lon,pole_lat=ecl_lonlat(polemean)
    pole_ra,pole_dec=ecl_to_equatorial(polemean)

    varpi_mean,varpi_R=circular_mean_deg([o["varpi_deg"] for o in sub])
    scatter=[angle_deg(v,pmean) for v in pvec]

    loo=leave_one_out(pvec,names,anti)
    shifts=[r["axis_shift_deg"] for r in loo]

    # Relation graph summary: nearest-neighbour / pairwise angular separations.
    pair=[]
    for i in range(len(sub)):
        for j in range(i+1,len(sub)):
            pair.append({
                "a":names[i],"b":names[j],
                "perihelion_angle_deg":angle_deg(pvec[i],pvec[j]),
                "pole_angle_deg":angle_deg(poles[i],poles[j])
            })
    pair_angles=[x["perihelion_angle_deg"] for x in pair]

    return {
        "status":"PASS",
        "n":len(sub),
        "cuts":spec,
        "objects":[{
            "name":o["name"],"a":o["a"],"q":o["q"],"e":o["e"],"i":o["i"],
            "om":o["om"],"w":o["w"],"varpi_deg":o["varpi_deg"],
            "condition_code":o["condition_code"],"data_arc_days":o["data_arc_days"],"n_obs":o["n_obs"]
        } for o in sub],
        "observed_relations":{
            "mean_perihelion_vector_ecl":pmean,
            "mean_perihelion_ecl_lon_deg":peri_lon,
            "mean_perihelion_ecl_lat_deg":peri_lat,
            "perihelion_vector_resultant":pR,
            "perihelion_scatter_deg_median":statistics.median(scatter),
            "perihelion_scatter_deg_p68":percentile(scatter,0.68),
            "circular_mean_varpi_deg":varpi_mean,
            "circular_varpi_resultant":varpi_R,
            "mean_orbital_pole_ecl":polemean,
            "mean_orbital_pole_ecl_lon_deg":pole_lon,
            "mean_orbital_pole_ecl_lat_deg":pole_lat,
            "mean_orbital_pole_ra_deg":pole_ra,
            "mean_orbital_pole_dec_deg":pole_dec,
            "pole_resultant":poleR,
            "pairwise_perihelion_angle_deg_median":statistics.median(pair_angles),
        },
        "relational_anti_apsidal_axis":{
            "interpretation":"geometric anti-axis of the observed perihelion-vector cluster; not instantaneous PBH position",
            "ecl_lon_deg":anti_lon,
            "ecl_lat_deg":anti_lat,
            "ra_deg":anti_ra,
            "dec_deg":anti_dec,
            "bootstrap":bootstrap_axis(pvec,anti),
            "leave_one_out_axis_shift_deg_max":max(shifts) if shifts else None,
            "leave_one_out_axis_shift_deg_median":statistics.median(shifts) if shifts else None,
        },
        "pairwise_relation_graph":pair,
    }

def main():
    signature,url,rows=fetch_rows()
    analyses={name:analyze_sample(rows,spec) for name,spec in SAMPLE_DEFS.items()}
    summary={
        "schema":"PBH_RELATIONAL_ORBIT_AXIS_V0_1",
        "status":"ANALYZED",
        "source":signature,
        "source_url":url,
        "base_query_count":len(rows),
        "samples":{},
        "method_boundary":[
            "JPL osculating elements are observations/model fits, not direct force residual vectors",
            "anti-apsidal axis is a relational geometric descriptor",
            "instantaneous PBH sky position, mass, distance and orbital phase are NOT identified by this step",
            "survey-selection bias is not corrected in v0.1"
        ],
        "canon_allowed":False
    }
    for name,a in analyses.items():
        if a.get("status")=="PASS":
            axis=a["relational_anti_apsidal_axis"]
            summary["samples"][name]={
                "n":a["n"],
                "anti_axis_ra_deg":axis["ra_deg"],
                "anti_axis_dec_deg":axis["dec_deg"],
                "anti_axis_ecl_lon_deg":axis["ecl_lon_deg"],
                "anti_axis_ecl_lat_deg":axis["ecl_lat_deg"],
                "bootstrap_p68_deg":axis["bootstrap"]["angular_offset_deg_p68"],
                "loo_shift_max_deg":axis["leave_one_out_axis_shift_deg_max"],
                "perihelion_resultant":a["observed_relations"]["perihelion_vector_resultant"],
                "varpi_mean_deg":a["observed_relations"]["circular_mean_varpi_deg"],
                "varpi_resultant":a["observed_relations"]["circular_varpi_resultant"],
                "objects":[o["name"] for o in a["objects"]]
            }
        else:
            summary["samples"][name]=a

    receipt={
        **summary,
        "full_analyses":analyses
    }
    print("PBH_ORBIT_AXIS_SUMMARY="+json.dumps(summary,sort_keys=True,separators=(",",":")),flush=True)
    print("PBH_ORBIT_AXIS_RECEIPT="+json.dumps(receipt,sort_keys=True,separators=(",",":")),flush=True)
    return 0

if __name__=="__main__":
    raise SystemExit(main())
