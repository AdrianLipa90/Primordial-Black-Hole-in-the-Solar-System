#!/usr/bin/env python3
from __future__ import annotations
import json, requests

API="https://data.minorplanetcenter.net/api/get-orb"
OBJECTS=["90377","2012 VP113","541132","2017 OF201","2013 SY99"]

out={"schema":"PBH_MPC_ORBIT_FIT_PROBE_V0_1","status":"ANALYZED","objects":{},"authority":"MPC_ORBIT_FIT_METADATA","canon_allowed":False}

for desig in OBJECTS:
    r=requests.get(API,json={"desig":desig},timeout=120,headers={"User-Agent":"GREMLIN-FERRET-PBH-ORBITFIT/0.1 read-only"})
    row={"http_status":r.status_code}
    if not r.ok:
        row["error"]=r.text[:500]; out["objects"][desig]=row; continue
    payload=r.json()
    if not payload or not payload[0].get("mpc_orb"):
        row["error"]="no mpc_orb"; out["objects"][desig]=row; continue
    orb=payload[0]["mpc_orb"][0]
    row["top_keys"]=sorted(orb.keys())
    for key in ["epoch_data","orbit_fit_statistics","designation_data","CAR","KEP","covariance","metadata"]:
        if key in orb: row[key]=orb[key]
    out["objects"][desig]=row

print("PBH_MPC_ORBIT_FIT_SUMMARY="+json.dumps(out["objects"],sort_keys=True,separators=(",",":")),flush=True)
print("PBH_MPC_ORBIT_FIT_RECEIPT="+json.dumps(out,sort_keys=True,separators=(",",":")),flush=True)
