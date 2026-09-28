#!/usr/bin/env python3
from __future__ import annotations

import json
import urllib.parse

import requests

API="https://ssd-api.jpl.nasa.gov/sbdb_query.api"
CONSTRAINT={
  "AND":[
    "a|GE|150",
    "q|GE|40",
    "e|LT|1",
    "data_arc|GE|365",
    "n_obs_used|GE|5"
  ]
}
FIELDS=[
  "spkid","pdes","full_name","class","orbit_id","epoch","e","a","q","i","om","w","ma",
  "sigma_e","sigma_a","sigma_q","sigma_i","sigma_om","sigma_w","sigma_ma",
  "data_arc","n_obs_used","condition_code","source"
]

params={
  "fields":",".join(FIELDS),
  "sb-kind":"a",
  "sb-class":"TNO",
  "sb-cdata":json.dumps(CONSTRAINT,separators=(",",":")),
  "full-prec":"true",
  "sort":"a"
}
url=API+"?"+urllib.parse.urlencode(params)
r=requests.get(url,timeout=120,headers={"User-Agent":"GREMLIN-FERRET-PBH-ORBITAL/0.1 read-only"})
r.raise_for_status()
data=r.json()

out={
  "schema":"PBH_JPL_EXTREME_TNO_INVENTORY_V0_1",
  "status":"ANALYZED",
  "source":"JPL_SBDB_QUERY_API",
  "source_url":url,
  "selection":{
    "a_min_au":150,
    "q_min_au":40,
    "e_max_exclusive":1,
    "data_arc_min_days":365,
    "n_obs_used_min":5,
    "class":"TNO"
  },
  "fields":data.get("fields"),
  "count":int(data.get("count",0)),
  "data":data.get("data",[]),
  "authority":"OBSERVATIONAL_ORBIT_CATALOG_ONLY",
  "canon_allowed":False
}
print("PBH_TNO_INVENTORY_SUMMARY="+json.dumps({
  "count":out["count"],
  "selection":out["selection"],
  "first_rows":out["data"][:12],
  "last_rows":out["data"][-5:],
},separators=(",",":"),sort_keys=True),flush=True)
print("PBH_TNO_INVENTORY_RECEIPT="+json.dumps(out,separators=(",",":"),sort_keys=True),flush=True)
