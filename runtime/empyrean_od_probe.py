#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
import requests

DESIG="2012 VP113"
OBS_API="https://data.minorplanetcenter.net/api/get-obs"

out={"schema":"PBH_EMPYREAN_OD_PROBE_V0_1","status":"STARTED","object":DESIG}

import empyrean

out["empyrean_version"]=getattr(empyrean,"version_string",lambda:"unknown")()
print("OD_PROBE_VERSION="+str(out["empyrean_version"]),flush=True)

# Initializes public SPICE/obscode/debiasing runtime data.
empyrean.initialize()
print("OD_PROBE_INIT=PASS",flush=True)

r=requests.get(
    OBS_API,
    json={"desigs":[DESIG],"output_format":["XML"],"ades_version":"2022"},
    timeout=120,
    headers={"User-Agent":"GREMLIN-FERRET-PBH-OD/0.1 read-only"},
)
r.raise_for_status()
payload=r.json()[0]["XML"]
path=Path("/tmp/2012_vp113.xml")
path.write_text(payload,encoding="utf-8")
out["mpc_xml_bytes"]=path.stat().st_size

obs,radar=empyrean.read_ades(str(path))
out["n_obs"]=len(obs)
out["radar_present"]=radar is not None and len(radar)>0
out["obs_type"]=type(obs).__name__
out["obs_attrs"]=[x for x in dir(obs) if not x.startswith("_")][:200]

for attr in ("stn","coordinates","obs_time","obstime","ra","dec","perm_id","prov_id","trk_sub"):
    try:
        value=getattr(obs,attr)
        row={"type":type(value).__name__}
        if hasattr(value,"to_pylist"):
            row["sample"]=value.to_pylist()[:5]
        elif hasattr(value,"to_numpy"):
            arr=value.to_numpy()
            row["sample"]=arr[:5].tolist()
        else:
            row["repr"]=repr(value)[:800]
        out.setdefault("obs_fields",{})[attr]=row
    except Exception as exc:
        out.setdefault("obs_fields",{})[attr]={"error":f"{type(exc).__name__}:{exc}"}

print("OD_PROBE_READ="+json.dumps({
    "n_obs":out["n_obs"],
    "radar_present":out["radar_present"],
    "obs_type":out["obs_type"],
    "obs_fields":out.get("obs_fields"),
},sort_keys=True,separators=(",",":")),flush=True)

fits=empyrean.determine(obs,radar=radar)
out["n_fits"]=len(fits)
out["n_delivered"]=len(fits.delivered)
out["failure_keys"]=list(fits.failures.keys())

def table_sample(x,n=8):
    if hasattr(x,"to_dataframe"):
        df=x.to_dataframe()
        return {"columns":[str(c) for c in df.columns],"rows":df.head(n).to_dict(orient="records")}
    if hasattr(x,"to_pyarrow"):
        t=x.to_pyarrow()
        return {"columns":t.column_names,"rows":t.slice(0,n).to_pylist()}
    return {"repr":repr(x)[:4000]}

out["summary"]=table_sample(fits.summary,5)
out["residuals"]=table_sample(fits.residuals,12)
out["orbits"]=table_sample(fits.orbits,3)
out["status"]="ANALYZED"

print("PBH_EMPYREAN_OD_PROBE_SUMMARY="+json.dumps({
    "version":out["empyrean_version"],
    "n_obs":out["n_obs"],
    "n_fits":out["n_fits"],
    "n_delivered":out["n_delivered"],
    "failures":out["failure_keys"],
    "summary":out["summary"],
    "residuals":out["residuals"],
    "orbits":out["orbits"],
},sort_keys=True,separators=(",",":"),default=str),flush=True)
print("PBH_EMPYREAN_OD_PROBE_RECEIPT="+json.dumps(out,sort_keys=True,separators=(",",":"),default=str),flush=True)
