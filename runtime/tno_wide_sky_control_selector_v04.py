#!/usr/bin/env python3
from __future__ import annotations

import json
import math
import urllib.parse

import requests

JPL_QUERY="https://ssd-api.jpl.nasa.gov/sbdb_query.api"
FREEZE_EPOCH_JD=2461311.5
N_SELECT=8
EXCLUDED_PDES={"90377","2012 VP113","541132"}
FIELDS=[
    "pdes","spkid","full_name","epoch","a","e","q","i","om","w","ma","n",
    "condition_code","data_arc","n_obs_used","two_body","source","last_obs"
]


def fnum(v):
    if v in (None,"","null"):
        return None
    try:
        x=float(v)
    except Exception:
        return None
    return x if math.isfinite(x) else None


def inum(v):
    x=fnum(v)
    return None if x is None else int(x)


def unit(v):
    s=math.sqrt(sum(x*x for x in v))
    if s<=0:
        raise ValueError("zero vector")
    return tuple(x/s for x in v)


def dot(a,b):
    return sum(x*y for x,y in zip(a,b))


def angle_deg(a,b):
    return math.degrees(math.acos(max(-1.0,min(1.0,dot(a,b)))))


def solve_kepler(M,e):
    E=M if e<0.8 else math.pi
    for _ in range(80):
        f=E-e*math.sin(E)-M
        fp=1-e*math.cos(E)
        d=f/fp
        E-=d
        if abs(d)<1e-14:
            break
    return E


def direction_from_elements(o):
    a=o["a"]; e=o["e"]; inc=math.radians(o["i"])
    Om=math.radians(o["om"]); w=math.radians(o["w"])
    M0=math.radians(o["ma"])
    nrad=math.radians(o["n"])
    M=(M0+nrad*(FREEZE_EPOCH_JD-o["epoch"]))%(2*math.pi)
    E=solve_kepler(M,e)
    cf=(math.cos(E)-e)/(1-e*math.cos(E))
    sf=(math.sqrt(max(0.0,1-e*e))*math.sin(E))/(1-e*math.cos(E))
    f=math.atan2(sf,cf)
    r=a*(1-e*math.cos(E))
    u=w+f
    x=r*(math.cos(Om)*math.cos(u)-math.sin(Om)*math.sin(u)*math.cos(inc))
    y=r*(math.sin(Om)*math.cos(u)+math.cos(Om)*math.sin(u)*math.cos(inc))
    z=r*(math.sin(u)*math.sin(inc))
    return unit((x,y,z))


def ecl_to_equatorial(v):
    eps=math.radians(23.439291111)
    x,y,z=v
    q=unit((x,y*math.cos(eps)-z*math.sin(eps),y*math.sin(eps)+z*math.cos(eps)))
    ra=math.degrees(math.atan2(q[1],q[0]))%360.0
    dec=math.degrees(math.asin(max(-1.0,min(1.0,q[2]))))
    return q,ra,dec


def fetch_pool():
    cdata=json.dumps({"AND":["a|GT|30","q|GT|25"]},separators=(",",":"))
    params={
        "fields":",".join(FIELDS),
        "sb-class":"TNO",
        "sb-cdata":cdata,
        "full-prec":"1",
        "limit":"5000",
    }
    response=requests.get(JPL_QUERY,params=params,timeout=90)
    response.raise_for_status()
    payload=response.json()
    fields=payload["fields"]
    rows=[]
    for raw in payload.get("data",[]):
        d=dict(zip(fields,raw))
        pdes=str(d.get("pdes") or "").strip()
        if pdes in EXCLUDED_PDES:
            continue
        vals={k:fnum(d.get(k)) for k in ("epoch","a","e","q","i","om","w","ma","n")}
        if any(v is None for v in vals.values()):
            continue
        cc=inum(d.get("condition_code"))
        arc=fnum(d.get("data_arc"))
        nobs=inum(d.get("n_obs_used"))
        two=str(d.get("two_body") or "").strip().upper()
        if cc is None or cc>2:
            continue
        if arc is None or arc<3652.5:
            continue
        if nobs is None or nobs<100:
            continue
        if two in {"Y","1","TRUE"}:
            continue
        o={
            "pdes":pdes,
            "spkid":str(d.get("spkid") or "").strip(),
            "full_name":str(d.get("full_name") or "").strip(),
            **vals,
            "condition_code":cc,
            "data_arc_days":arc,
            "n_obs_used":nobs,
            "source":d.get("source"),
            "last_obs":d.get("last_obs"),
        }
        ecl=direction_from_elements(o)
        equ,ra,dec=ecl_to_equatorial(ecl)
        o["approx_ecl_unit"]=list(ecl)
        o["approx_equ_unit"]=list(equ)
        o["approx_ra_deg"]=ra
        o["approx_dec_deg"]=dec
        rows.append(o)
    rows.sort(key=lambda o:o["pdes"])
    return payload.get("signature"),response.url,rows


def select_wide_sky(pool,n):
    if len(pool)<n:
        raise ValueError("eligible pool smaller than requested control sample")

    best=None
    for i in range(len(pool)):
        for j in range(i+1,len(pool)):
            sep=angle_deg(pool[i]["approx_equ_unit"],pool[j]["approx_equ_unit"])
            key=(sep,pool[i]["pdes"],pool[j]["pdes"])
            if best is None or key[0]>best[0] or (abs(key[0]-best[0])<1e-12 and key[1:]<best[1:]):
                best=(sep,pool[i]["pdes"],pool[j]["pdes"],i,j)

    selected=[best[3],best[4]]
    trace=[{
        "step":1,
        "pdes":[pool[best[3]]["pdes"],pool[best[4]]["pdes"]],
        "pair_separation_deg":best[0],
    }]

    while len(selected)<n:
        choice=None
        for i,o in enumerate(pool):
            if i in selected:
                continue
            minsep=min(angle_deg(o["approx_equ_unit"],pool[j]["approx_equ_unit"]) for j in selected)
            meansep=sum(angle_deg(o["approx_equ_unit"],pool[j]["approx_equ_unit"]) for j in selected)/len(selected)
            key=(minsep,meansep,o["pdes"])
            if choice is None or key[0]>choice[0] or (
                abs(key[0]-choice[0])<1e-12 and (
                    key[1]>choice[1] or (
                        abs(key[1]-choice[1])<1e-12 and key[2]<choice[2]
                    )
                )
            ):
                choice=(minsep,meansep,o["pdes"],i)
        selected.append(choice[3])
        trace.append({
            "step":len(selected)-1,
            "pdes":pool[choice[3]]["pdes"],
            "min_separation_deg_when_added":choice[0],
            "mean_separation_deg_when_added":choice[1],
        })

    return [pool[i] for i in selected],trace


def main():
    signature,url,pool=fetch_pool()
    selected,trace=select_wide_sky(pool,N_SELECT)
    pairwise=[]
    for i in range(len(selected)):
        for j in range(i+1,len(selected)):
            pairwise.append(angle_deg(selected[i]["approx_equ_unit"],selected[j]["approx_equ_unit"]))

    receipt={
        "schema":"PBH_WIDE_SKY_TNO_CONTROL_SELECTION_V0_4",
        "status":"FROZEN_BEFORE_RESIDUALS",
        "authority":"JPL_SBDB_SELECTION_ONLY",
        "freeze_epoch_jd":FREEZE_EPOCH_JD,
        "source_signature":signature,
        "query_url":url,
        "criteria":{
            "sb_class":"TNO",
            "a_gt_au":30.0,
            "q_gt_au":25.0,
            "data_arc_ge_days":3652.5,
            "n_obs_used_ge":100,
            "condition_code_le":2,
            "two_body_rejected":True,
            "excluded_target_pdes":sorted(EXCLUDED_PDES),
            "n_controls":N_SELECT,
        },
        "eligible_count":len(pool),
        "selection_algorithm":"max-separation pair, then deterministic farthest-point max-min on approximate propagated J2000 sky directions",
        "selection_trace":trace,
        "selected":selected,
        "pairwise_separation_deg":{
            "min":min(pairwise),
            "median":sorted(pairwise)[len(pairwise)//2],
            "max":max(pairwise),
        },
        "anti_confirmation":{
            "residuals_seen_before_selection":False,
            "historical_PBH_coordinate_used":False,
            "three_tracer_STF_axis_used":False,
            "historical_576_AU_used":False,
        },
        "canon_allowed":False,
    }
    print("PBH_TNO_CONTROL_SELECTION="+json.dumps(receipt,sort_keys=True),flush=True)


if __name__=="__main__":
    main()
