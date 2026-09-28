#!/usr/bin/env python3
from __future__ import annotations
import json, requests

API="https://data.minorplanetcenter.net/api/get-obs"
OBJECTS=["90377","2012 VP113","541132","2017 OF201","2013 SY99"]

out={"schema":"PBH_MPC_ASTROMETRY_PROBE_V0_1","status":"ANALYZED","objects":{},"authority":"RAW_ASTROMETRY_ONLY","canon_allowed":False}
for desig in OBJECTS:
    r=requests.get(API,json={"desigs":[desig],"output_format":["ADES_DF"],"ades_version":"2022"},timeout=120,
                   headers={"User-Agent":"GREMLIN-FERRET-PBH-ASTROMETRY/0.1 read-only"})
    row={"http_status":r.status_code}
    if not r.ok:
        row["error"]=r.text[:500]
        out["objects"][desig]=row
        continue
    data=r.json()[0]["ADES_DF"]
    row["count"]=len(data)
    row["columns"]=sorted({k for rec in data for k in rec.keys()})
    row["first"]=data[0] if data else None
    row["last"]=data[-1] if data else None
    out["objects"][desig]=row

print("PBH_MPC_ASTROMETRY_PROBE_SUMMARY="+json.dumps({
  k:{"count":v.get("count"),"columns":v.get("columns"),"first":v.get("first"),"last":v.get("last"),"error":v.get("error")}
  for k,v in out["objects"].items()
},sort_keys=True,separators=(",",":")),flush=True)
print("PBH_MPC_ASTROMETRY_PROBE_RECEIPT="+json.dumps(out,sort_keys=True,separators=(",",":")),flush=True)
