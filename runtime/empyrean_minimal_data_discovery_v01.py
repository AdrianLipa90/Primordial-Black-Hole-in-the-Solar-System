#!/usr/bin/env python3
from __future__ import annotations
import json, os, site
from pathlib import Path

def main():
    roots=[Path(p) for p in site.getsitepackages()]
    roots.append(Path(site.getusersitepackages()))
    pats=("de440","gm","naif")
    hits=[]
    for root in roots:
        if not root.exists():
            continue
        for p in root.rglob("*"):
            if not p.is_file():
                continue
            name=p.name.lower()
            if (
                "de440" in name
                or name in {"gm_de440.tpc","gm_de431.tpc","gm.tpc"}
                or (name.startswith("gm") and p.suffix.lower() in {".tpc",".txt",".json"})
            ):
                hits.append({"path":str(p),"size":p.stat().st_size})
    print("PBH_EMPYREAN_DATA_FILES="+json.dumps({
        "schema":"PBH_EMPYREAN_MINIMAL_DATA_DISCOVERY_V0_1",
        "hits":sorted(hits,key=lambda x:x["path"]),
    },sort_keys=True), flush=True)

if __name__=="__main__":
    main()
