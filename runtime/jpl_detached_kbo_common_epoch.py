#!/usr/bin/env python3
from __future__ import annotations
import json, urllib.parse, urllib.request

OBJECTS=[
"2000 CR105","Sedna","2004 VN112","2010 GB174","2012 VP113",
"2013 FT28","2013 RA109","2013 SY99","2013 UT15","2014 SR349","2015 RX245"
]

def fetch(name):
    url="https://ssd-api.jpl.nasa.gov/sbdb.api?"+urllib.parse.urlencode(
        {"sstr":name,"full-prec":"1","cd-epoch":"1"}
    )
    req=urllib.request.Request(url,headers={"User-Agent":"GREMLIN-PBH-common-epoch/0.1 read-only"})
    with urllib.request.urlopen(req,timeout=30) as resp:
        d=json.loads(resp.read())
    o=d["orbit"]
    elems={e["name"]:e for e in o["elements"]}
    out={
        "query":name,"fullname":d["object"].get("fullname"),"orbit_id":o.get("orbit_id"),
        "epoch":float(o["epoch"]),"cd_epoch":o.get("cd_epoch"),"equinox":o.get("equinox"),
        "source_url":url,
    }
    for key in ("a","e","q","i","om","w","ma","n","per","ad"):
        if key in elems:
            out[key]=float(elems[key]["value"])
            out[key+"_sigma"]=None if elems[key].get("sigma") is None else float(elems[key]["sigma"])
    return out

def main():
    rows=[fetch(x) for x in OBJECTS]
    epochs=sorted(set(r["epoch"] for r in rows))
    receipt={
      "schema":"PBH_JPL_DETACHED_KBO_COMMON_EPOCH_V0_1",
      "status":"PASS" if len(epochs)==1 else "FAIL_EPOCH_MISMATCH",
      "epochs":epochs,
      "objects":rows,
      "authority":"DATA_RECOVERY_ONLY",
      "canon_allowed":False,
    }
    print("PBH_COMMON_EPOCH="+json.dumps(receipt,sort_keys=True,separators=(",",":")),flush=True)
    return 0 if len(epochs)==1 else 2

if __name__=="__main__": raise SystemExit(main())
