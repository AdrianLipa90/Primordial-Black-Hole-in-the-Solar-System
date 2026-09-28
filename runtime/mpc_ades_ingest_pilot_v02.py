#!/usr/bin/env python3
from __future__ import annotations

import collections
import hashlib
import json
import math
from datetime import datetime

import requests

MPC_URL = "https://data.minorplanetcenter.net/api/get-obs"
OBJECTS = ["Sedna", "2012 VP113", "2015 TG387"]


def numeric_summary(rows, key):
    vals = []
    for row in rows:
        try:
            v = float(row.get(key))
        except Exception:
            continue
        if math.isfinite(v):
            vals.append(v)
    if not vals:
        return None
    mean = sum(vals) / len(vals)
    var = sum((x - mean) ** 2 for x in vals) / max(len(vals) - 1, 1)
    return {
        "n": len(vals),
        "min": min(vals),
        "max": max(vals),
        "mean": mean,
        "std": math.sqrt(var),
    }


def counts(rows, key, top=20):
    c = collections.Counter(str(r.get(key, "")) for r in rows if r.get(key) not in (None, ""))
    return c.most_common(top)


def ingest_one(desig):
    response = requests.get(
        MPC_URL,
        json={"desigs": [desig], "output_format": ["ADES_DF"], "ades_version": "2022"},
        timeout=90,
    )
    response.raise_for_status()
    payload = response.json()
    if not payload or "ADES_DF" not in payload[0]:
        raise RuntimeError(f"{desig}: ADES_DF missing")
    rows = payload[0]["ADES_DF"]
    canonical = json.dumps(rows, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode("utf-8")
    keys = sorted({str(k) for row in rows for k in row.keys()})

    times = sorted(str(r.get("obsTime")) for r in rows if r.get("obsTime"))
    receipt = {
        "designation": desig,
        "status": "INGESTED",
        "n_rows": len(rows),
        "sha256_canonical_ades_df": hashlib.sha256(canonical).hexdigest(),
        "keys": keys,
        "first_obsTime": times[0] if times else None,
        "last_obsTime": times[-1] if times else None,
        "top_stations": counts(rows, "stn"),
        "top_modes": counts(rows, "mode"),
        "top_catalogs": counts(rows, "astCat"),
        "sigRA": numeric_summary(rows, "sigRA"),
        "sigDec": numeric_summary(rows, "sigDec"),
        "rmsRA": numeric_summary(rows, "rmsRA"),
        "rmsDec": numeric_summary(rows, "rmsDec"),
        "ra": numeric_summary(rows, "ra"),
        "dec": numeric_summary(rows, "dec"),
        "sample_rows": rows[:2],
    }
    return receipt


def main():
    results = []
    failures = []
    for desig in OBJECTS:
        try:
            r = ingest_one(desig)
            results.append(r)
            print("PBH_ADES_OBJECT=" + json.dumps(r, ensure_ascii=False, sort_keys=True), flush=True)
        except Exception as exc:
            f = {"designation": desig, "error_type": type(exc).__name__, "error": str(exc)}
            failures.append(f)
            print("PBH_ADES_FAILURE=" + json.dumps(f, sort_keys=True), flush=True)

    final = {
        "schema": "PBH_MPC_ADES_INGEST_V0_2",
        "status": "PASS" if results and not failures else ("PARTIAL" if results else "FAIL"),
        "authority": "OBSERVATIONAL_DATA_ONLY",
        "mpc_endpoint": MPC_URL,
        "ades_version": "2022",
        "results": results,
        "failures": failures,
        "canon_allowed": False,
    }
    print("PBH_ADES_RECEIPT=" + json.dumps(final, ensure_ascii=False, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
