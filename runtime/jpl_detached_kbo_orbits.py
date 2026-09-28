#!/usr/bin/env python3
from __future__ import annotations

import json
import urllib.parse
import urllib.request

OBJECTS = [
    "2000 CR105",
    "90377 Sedna",
    "2004 VN112",
    "2010 GB174",
    "2012 VP113",
    "2013 FT28",
    "2013 RA109",
    "2013 SY99",
    "2013 UT15",
    "2014 SR349",
    "2015 RX245",
]


def fetch(name: str) -> dict:
    params = urllib.parse.urlencode(
        {
            "sstr": name,
            "full-prec": "1",
            "cov": "mat",
            "soln-epoch": "1",
            "cd-epoch": "1",
        }
    )
    url = "https://ssd-api.jpl.nasa.gov/sbdb.api?" + params
    req = urllib.request.Request(
        url,
        headers={"User-Agent": "GREMLIN-PBH-orbital-recovery/0.1 read-only"},
    )
    with urllib.request.urlopen(req, timeout=30) as resp:
        payload = resp.read()
    data = json.loads(payload)

    orbit = data["orbit"]
    elements = {
        e["name"]: {
            "value": e.get("value"),
            "sigma": e.get("sigma"),
            "units": e.get("units"),
        }
        for e in orbit.get("elements", [])
    }

    covariance = orbit.get("covariance") or {}
    return {
        "query": name,
        "source_url": url,
        "fullname": data.get("object", {}).get("fullname"),
        "spkid": data.get("object", {}).get("spkid"),
        "orbit_id": orbit.get("orbit_id"),
        "epoch": orbit.get("epoch"),
        "cd_epoch": orbit.get("cd_epoch"),
        "equinox": orbit.get("equinox"),
        "rms": orbit.get("rms"),
        "condition_code": orbit.get("condition_code"),
        "data_arc": orbit.get("data_arc"),
        "n_obs_used": orbit.get("n_obs_used"),
        "first_obs": orbit.get("first_obs"),
        "last_obs": orbit.get("last_obs"),
        "source": orbit.get("source"),
        "producer": orbit.get("producer"),
        "elements": elements,
        "covariance_labels": covariance.get("labels"),
        "covariance": covariance.get("data"),
        "cov_epoch": covariance.get("epoch"),
    }


def main() -> int:
    rows = []
    for name in OBJECTS:
        try:
            row = fetch(name)
            print(
                "JPL_FETCH_OK="
                + json.dumps(
                    {
                        "query": name,
                        "fullname": row.get("fullname"),
                        "orbit_id": row.get("orbit_id"),
                        "n_obs_used": row.get("n_obs_used"),
                    },
                    sort_keys=True,
                    separators=(",", ":"),
                ),
                flush=True,
            )
            rows.append(row)
        except Exception as exc:
            rows.append({"query": name, "error": f"{type(exc).__name__}: {exc}"})
            print(
                "JPL_FETCH_FAIL="
                + json.dumps({"query": name, "error": rows[-1]["error"]}),
                flush=True,
            )

    receipt = {
        "schema": "PBH_JPL_DETACHED_KBO_ORBITS_V0_1",
        "status": "ANALYZED"
        if all("error" not in row for row in rows)
        else "PARTIAL",
        "selection": {
            "origin": "Brown_Batygin_2021_Table_1",
            "criteria": "150<a<1000 AU and q>42 AU, multi-opposition, as of 2021-08-20",
            "objects": OBJECTS,
        },
        "source": {
            "authority": "JPL_SSD_SBDB_API",
            "api": "https://ssd-api.jpl.nasa.gov/sbdb.api",
            "full_precision": True,
            "covariance": "mat",
            "solution_epoch": True,
        },
        "objects": rows,
        "interpretation": "ORBIT_DATA_RECOVERY_ONLY",
        "canon_allowed": False,
    }
    print(
        "PBH_JPL_ORBITS_RECEIPT="
        + json.dumps(receipt, sort_keys=True, separators=(",", ":")),
        flush=True,
    )
    return 0 if receipt["status"] == "ANALYZED" else 2


if __name__ == "__main__":
    raise SystemExit(main())
