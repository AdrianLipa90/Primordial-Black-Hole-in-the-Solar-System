#!/usr/bin/env python3
import json, requests

URL="https://ssd.jpl.nasa.gov/api/horizons.api"
tests=[
    {
        "label":"sedna_261_1990",
        "params":{
            "format":"json",
            "COMMAND":"'90377'",
            "OBJ_DATA":"'NO'",
            "MAKE_EPHEM":"'YES'",
            "EPHEM_TYPE":"'OBSERVER'",
            "CENTER":"'261'",
            "TLIST":"'1990-09-25 10:34:30.432'",
            "TIME_TYPE":"'UT'",
            "TIME_DIGITS":"'FRACSEC'",
            "QUANTITIES":"'1'",
            "CSV_FORMAT":"'YES'",
        },
    },
    {
        "label":"sedna_568_2004",
        "params":{
            "format":"json",
            "COMMAND":"'90377'",
            "OBJ_DATA":"'NO'",
            "MAKE_EPHEM":"'YES'",
            "EPHEM_TYPE":"'OBSERVER'",
            "CENTER":"'568'",
            "TLIST":"'2004-01-01 00:00:00'",
            "TIME_TYPE":"'UT'",
            "TIME_DIGITS":"'FRACSEC'",
            "QUANTITIES":"'1'",
            "CSV_FORMAT":"'YES'",
        },
    },
]
for t in tests:
    try:
        r=requests.get(URL,params=t["params"],timeout=60)
        out={"label":t["label"],"http_status":r.status_code,"url":r.url}
        try:
            payload=r.json()
            out["signature"]=payload.get("signature")
            result=payload.get("result","")
            out["result_head"]=result[:2500]
            out["has_ephemeris_markers"]="$$SOE" in result and "$$EOE" in result
        except Exception:
            out["text_head"]=r.text[:2500]
        print("PBH_HORIZONS_CENTER_PROBE="+json.dumps(out,sort_keys=True),flush=True)
    except Exception as exc:
        print("PBH_HORIZONS_CENTER_PROBE="+json.dumps({"label":t["label"],"error_type":type(exc).__name__,"error":str(exc)},sort_keys=True),flush=True)
