#!/usr/bin/env python3
from __future__ import annotations

import csv
import hashlib
import io
import json
import math
import statistics
from collections import defaultdict
from datetime import datetime, timezone

import requests

MPC_URL = "https://data.minorplanetcenter.net/api/get-obs"
HORIZONS_URL = "https://ssd.jpl.nasa.gov/api/horizons.api"

OBJECTS = {
    "2012 VP113": "'DES=2012 VP113;'",
    "2015 TG387": "'541132;'",
}


def as_float(x):
    try:
        v = float(x)
        return v if math.isfinite(v) else None
    except Exception:
        return None


def parse_iso_utc(s: str) -> datetime:
    return datetime.fromisoformat(s.replace("Z", "+00:00")).astimezone(timezone.utc)


def horizons_time(s: str) -> str:
    dt = parse_iso_utc(s)
    # Preserve millisecond-level timestamp used by MPC/Horizons.
    return dt.strftime("%Y-%b-%d %H:%M:%S.") + f"{dt.microsecond // 1000:03d}"


def fetch_mpc(desig: str):
    r = requests.get(
        MPC_URL,
        json={"desigs": [desig], "output_format": ["ADES_DF"], "ades_version": "2022"},
        timeout=90,
    )
    r.raise_for_status()
    payload = r.json()
    if not payload or "ADES_DF" not in payload[0]:
        raise RuntimeError(f"{desig}: MPC ADES_DF missing")
    rows = payload[0]["ADES_DF"]
    out = []
    for row in rows:
        ra = as_float(row.get("ra"))
        dec = as_float(row.get("dec"))
        stn = str(row.get("stn") or "").strip()
        obstime = str(row.get("obstime") or "").strip()
        if ra is None or dec is None or not stn or not obstime:
            continue
        out.append({
            "obs_id": row.get("obsid"),
            "time": obstime,
            "stn": stn,
            "catalog": row.get("astcat"),
            "mode": row.get("mode"),
            "ra_obs_deg": ra,
            "dec_obs_deg": dec,
            "rmsra_reported_arcsec": as_float(row.get("rmsra")),
            "rmsdec_reported_arcsec": as_float(row.get("rmsdec")),
        })
    return out


def parse_horizons_result(result: str):
    if "$$SOE" not in result or "$$EOE" not in result:
        raise RuntimeError("Horizons result has no SOE/EOE block")
    body = result.split("$$SOE", 1)[1].split("$$EOE", 1)[0].strip()
    rows = []
    for rec in csv.reader(io.StringIO(body)):
        vals = [x.strip() for x in rec]
        if len(vals) < 6:
            continue
        cal = vals[0]
        jd = as_float(vals[1])
        ra = as_float(vals[-2])
        dec = as_float(vals[-1])
        if jd is None or ra is None or dec is None:
            raise RuntimeError(f"unparseable Horizons row: {vals!r}")
        rows.append({"cal": cal, "jd_ut": jd, "ra_pred_deg": ra, "dec_pred_deg": dec})
    return rows


def query_horizons(command: str, stn: str, obs_rows):
    if not obs_rows:
        return []
    # Horizons accepts a quoted space-separated calendar TLIST.
    tlist = " ".join(f"'{horizons_time(x['time'])}'" for x in obs_rows)
    params = {
        "format": "json",
        "COMMAND": command,
        "OBJ_DATA": "'NO'",
        "MAKE_EPHEM": "'YES'",
        "EPHEM_TYPE": "'OBSERVER'",
        "CENTER": f"'{stn}@399'",
        "TLIST": tlist,
        "TLIST_TYPE": "'CAL'",
        "TIME_TYPE": "'UT'",
        "QUANTITIES": "'1'",
        "REF_SYSTEM": "'ICRF'",
        "ANG_FORMAT": "'DEG'",
        "CSV_FORMAT": "'YES'",
        "EXTRA_PREC": "'YES'",
        "CAL_FORMAT": "'BOTH'",
    }
    r = requests.get(HORIZONS_URL, params=params, timeout=90)
    r.raise_for_status()
    payload = r.json()
    result = payload.get("result", "")
    if "** Error" in result:
        raise RuntimeError(f"Horizons error for {stn}: {result[:500]}")
    pred = parse_horizons_result(result)
    if len(pred) != len(obs_rows):
        raise RuntimeError(
            f"{stn}: Horizons returned {len(pred)} rows for {len(obs_rows)} observations"
        )
    return pred


def wrap_deg(x):
    return ((x + 180.0) % 360.0) - 180.0


def residual_row(obj, obs, pred):
    dra = wrap_deg(obs["ra_obs_deg"] - pred["ra_pred_deg"])
    dra_as = dra * 3600.0 * math.cos(math.radians(pred["dec_pred_deg"]))
    ddec_as = (obs["dec_obs_deg"] - pred["dec_pred_deg"]) * 3600.0
    return {
        "object": obj,
        "obs_id": obs["obs_id"],
        "time": obs["time"],
        "jd_ut": pred["jd_ut"],
        "stn": obs["stn"],
        "catalog": obs["catalog"],
        "mode": obs["mode"],
        "ra_obs_deg": obs["ra_obs_deg"],
        "dec_obs_deg": obs["dec_obs_deg"],
        "ra_pred_deg": pred["ra_pred_deg"],
        "dec_pred_deg": pred["dec_pred_deg"],
        "dra_cosdec_arcsec": dra_as,
        "ddec_arcsec": ddec_as,
        "rho_arcsec": math.hypot(dra_as, ddec_as),
        "rmsra_reported_arcsec": obs["rmsra_reported_arcsec"],
        "rmsdec_reported_arcsec": obs["rmsdec_reported_arcsec"],
    }


def med_mad(vals):
    vals = [float(x) for x in vals if x is not None and math.isfinite(float(x))]
    if not vals:
        return None
    med = statistics.median(vals)
    mad = statistics.median(abs(x - med) for x in vals)
    return {
        "n": len(vals),
        "mean": sum(vals) / len(vals),
        "median": med,
        "mad": mad,
        "rms": math.sqrt(sum(x*x for x in vals) / len(vals)),
        "min": min(vals),
        "max": max(vals),
    }


def summarize(rows):
    return {
        "dra_cosdec_arcsec": med_mad([r["dra_cosdec_arcsec"] for r in rows]),
        "ddec_arcsec": med_mad([r["ddec_arcsec"] for r in rows]),
        "rho_arcsec": med_mad([r["rho_arcsec"] for r in rows]),
    }


def group_summary(rows, key, minimum=2):
    groups = defaultdict(list)
    for r in rows:
        groups[str(r.get(key) or "UNKNOWN")].append(r)
    out = {}
    for k, rr in groups.items():
        if len(rr) >= minimum:
            out[k] = {"n": len(rr), **summarize(rr)}
    return out


def solve3(a, b):
    m = [list(map(float, row)) + [float(b[i])] for i, row in enumerate(a)]
    for col in range(3):
        piv = max(range(col, 3), key=lambda i: abs(m[i][col]))
        if abs(m[piv][col]) < 1e-20:
            raise RuntimeError("singular quadratic normal matrix")
        m[col], m[piv] = m[piv], m[col]
        q = m[col][col]
        for j in range(col, 4):
            m[col][j] /= q
        for i in range(3):
            if i == col:
                continue
            q = m[i][col]
            for j in range(col, 4):
                m[i][j] -= q * m[col][j]
    return [m[i][3] for i in range(3)]


def quadratic_fit(rows, field):
    if len(rows) < 5:
        return None
    jds = [r["jd_ut"] for r in rows]
    t0 = sum(jds) / len(jds)
    ts = [(x - t0) / 365.25 for x in jds]
    ys = [r[field] for r in rows]
    s0 = len(ts)
    s1 = sum(ts)
    s2 = sum(t*t for t in ts)
    s3 = sum(t**3 for t in ts)
    s4 = sum(t**4 for t in ts)
    b0 = sum(ys)
    b1 = sum(t*y for t, y in zip(ts, ys))
    b2 = sum((t*t)*y for t, y in zip(ts, ys))
    c0, c1, c2 = solve3(
        [[s0, s1, s2], [s1, s2, s3], [s2, s3, s4]],
        [b0, b1, b2],
    )
    resid = [y - (c0 + c1*t + c2*t*t) for t, y in zip(ts, ys)]
    return {
        "t0_jd_ut": t0,
        "intercept_arcsec": c0,
        "linear_arcsec_per_year": c1,
        "quadratic_arcsec_per_year2": c2,
        "second_derivative_arcsec_per_year2": 2.0*c2,
        "postfit_rms_arcsec": math.sqrt(sum(x*x for x in resid) / len(resid)),
    }


def main():
    all_rows = []
    object_receipts = {}
    failures = []

    for obj, command in OBJECTS.items():
        try:
            obs = fetch_mpc(obj)
            by_station = defaultdict(list)
            for row in obs:
                by_station[row["stn"]].append(row)
            object_rows = []
            station_failures = []
            for stn in sorted(by_station):
                rr = sorted(by_station[stn], key=lambda x: x["time"])
                try:
                    pp = query_horizons(command, stn, rr)
                    for o, p in zip(rr, pp, strict=True):
                        object_rows.append(residual_row(obj, o, p))
                except Exception as exc:
                    station_failures.append({
                        "stn": stn,
                        "n": len(rr),
                        "error_type": type(exc).__name__,
                        "error": str(exc),
                    })
            object_rows.sort(key=lambda x: (x["jd_ut"], x["stn"], str(x["obs_id"])))
            all_rows.extend(object_rows)
            object_receipts[obj] = {
                "n_mpc_usable": len(obs),
                "n_residuals": len(object_rows),
                "n_station_failures": len(station_failures),
                "station_failures": station_failures,
                "summary": summarize(object_rows),
                "by_station": group_summary(object_rows, "stn", minimum=2),
                "by_catalog": group_summary(object_rows, "catalog", minimum=2),
                "quadratic_profile": {
                    "dra": quadratic_fit(object_rows, "dra_cosdec_arcsec"),
                    "ddec": quadratic_fit(object_rows, "ddec_arcsec"),
                },
            }
        except Exception as exc:
            failures.append({
                "object": obj,
                "error_type": type(exc).__name__,
                "error": str(exc),
            })

    canonical = json.dumps(all_rows, sort_keys=True, separators=(",", ":")).encode()
    receipt = {
        "schema": "PBH_MPC_HORIZONS_OC_RESIDUALS_V0_1",
        "status": (
            "PASS" if object_receipts and not failures
            and all(v["n_station_failures"] == 0 for v in object_receipts.values())
            else ("PARTIAL" if object_receipts else "FAIL")
        ),
        "authority": "OBSERVATIONAL_OC_ONLY",
        "mpc_endpoint": MPC_URL,
        "horizons_endpoint": HORIZONS_URL,
        "horizons_baseline": "current JPL solution; not statistically independent of MPC astrometry",
        "objects": object_receipts,
        "n_rows": len(all_rows),
        "rows_sha256": hashlib.sha256(canonical).hexdigest(),
        "failures": failures,
        "interpretation_boundary": {
            "not_a_PBH_detection": True,
            "not_independent_of_JPL_orbit_fit": True,
            "purpose": "validate exact-station exact-time O-C structure before a common-tide inversion",
        },
        "canon_allowed": False,
    }
    print("PBH_OC_RECEIPT=" + json.dumps(receipt, sort_keys=True), flush=True)
    for row in all_rows:
        print("PBH_OC_ROW=" + json.dumps(row, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
