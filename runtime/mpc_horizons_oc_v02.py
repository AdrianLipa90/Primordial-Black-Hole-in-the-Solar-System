#!/usr/bin/env python3
from __future__ import annotations

import collections
import hashlib
import json
import math
import statistics
import time
from datetime import datetime, timezone

import requests

MPC_URL = "https://data.minorplanetcenter.net/api/get-obs"
HORIZONS_URL = "https://ssd.jpl.nasa.gov/api/horizons.api"

OBJECTS = [
    {"designation": "Sedna", "command": "90377;"},
    {"designation": "2012 VP113", "command": "DES=2012 VP113;"},
    {"designation": "2015 TG387", "command": "541132;"},
]

CHUNK = 80


def ffloat(x):
    try:
        v = float(x)
    except Exception:
        return None
    return v if math.isfinite(v) else None


def iso_to_jd(s: str) -> float:
    dt = datetime.fromisoformat(s.replace("Z", "+00:00"))
    if dt.tzinfo is None:
        dt = dt.replace(tzinfo=timezone.utc)
    return dt.timestamp() / 86400.0 + 2440587.5


def wrap_deg(x: float) -> float:
    return (x + 180.0) % 360.0 - 180.0


def chunks(xs, n):
    for i in range(0, len(xs), n):
        yield xs[i:i+n]


def horizons_batch(command: str, station: str, rows: list[dict]) -> tuple[list[dict], dict]:
    jd_list = [iso_to_jd(r["obstime"]) for r in rows]
    tlist = ",".join(f"'{jd:.12f}'" for jd in jd_list)
    params = {
        "format": "json",
        "COMMAND": f"'{command}'",
        "OBJ_DATA": "'NO'",
        "MAKE_EPHEM": "'YES'",
        "EPHEM_TYPE": "'OBSERVER'",
        "CENTER": f"'{station}'",
        "TLIST": tlist,
        "TLIST_TYPE": "'JD'",
        "TIME_TYPE": "'UT'",
        "TIME_DIGITS": "'FRACSEC'",
        "CAL_FORMAT": "'JD'",
        "ANG_FORMAT": "'DEG'",
        "CSV_FORMAT": "'YES'",
        "QUANTITIES": "'1'",
        "SKIP_DAYLT": "'NO'",
        "ELEV_CUT": "'-90'",
        "APPARENT": "'AIRLESS'",
    }
    response = requests.get(HORIZONS_URL, params=params, timeout=90)
    response.raise_for_status()
    payload = response.json()
    result = payload.get("result", "")
    if "$$SOE" not in result or "$$EOE" not in result:
        raise RuntimeError("Horizons response lacks ephemeris markers: " + result[:500])

    block = result.split("$$SOE", 1)[1].split("$$EOE", 1)[0]
    pred = []
    for line in block.splitlines():
        line = line.strip()
        if not line:
            continue
        parts = [p.strip() for p in line.split(",")]
        nums = []
        for p in parts:
            try:
                nums.append(float(p))
            except Exception:
                pass
        if len(nums) < 3:
            continue
        pred.append({"jd": nums[0], "ra_deg": nums[1], "dec_deg": nums[2]})

    if len(pred) != len(rows):
        raise RuntimeError(
            f"Horizons row mismatch station={station}: expected {len(rows)}, got {len(pred)}; "
            f"data_block={block[:1000]}"
        )

    meta = {
        "signature": payload.get("signature"),
        "target_line": next((x.strip() for x in result.splitlines() if x.startswith("Target body name:")), None),
        "center_site_line": next((x.strip() for x in result.splitlines() if x.startswith("Center-site name:")), None),
        "initial_orbit_rms_line": next((x.strip() for x in result.splitlines() if "Residual RMS=" in x), None),
    }
    return pred, meta


def rms(vals):
    return math.sqrt(sum(x*x for x in vals) / len(vals)) if vals else None


def median(vals):
    return statistics.median(vals) if vals else None


def mad(vals):
    if not vals:
        return None
    m = statistics.median(vals)
    return statistics.median(abs(x-m) for x in vals)


def linear_slope_per_year(rows, key):
    pts = [(r["jd"], r[key]) for r in rows if math.isfinite(r[key])]
    if len(pts) < 3:
        return None
    t0 = sum(t for t, _ in pts) / len(pts)
    xs = [(t - t0) / 365.25 for t, _ in pts]
    ys = [y for _, y in pts]
    xm = sum(xs) / len(xs)
    ym = sum(ys) / len(ys)
    den = sum((x-xm)**2 for x in xs)
    if den <= 0:
        return None
    return sum((x-xm)*(y-ym) for x, y in zip(xs, ys)) / den


def demean_by_group(rows, key, group_fields=("stn", "astcat")):
    groups = collections.defaultdict(list)
    for idx, r in enumerate(rows):
        groups[tuple(r.get(k) or "" for k in group_fields)].append((idx, r[key]))
    out = [None] * len(rows)
    for vals in groups.values():
        mu = sum(v for _, v in vals) / len(vals)
        for idx, v in vals:
            out[idx] = v - mu
    return out


def summarize_residuals(rows):
    ra = [r["dra_cosdec_arcsec"] for r in rows]
    de = [r["ddec_arcsec"] for r in rows]
    ra_dm = demean_by_group(rows, "dra_cosdec_arcsec")
    de_dm = demean_by_group(rows, "ddec_arcsec")
    work = []
    for r, x, y in zip(rows, ra_dm, de_dm):
        q = dict(r)
        q["dra_cosdec_demeaned_arcsec"] = x
        q["ddec_demeaned_arcsec"] = y
        work.append(q)
    return {
        "n": len(rows),
        "dra_cosdec_arcsec": {
            "median": median(ra),
            "mad": mad(ra),
            "rms": rms(ra),
            "slope_arcsec_per_year": linear_slope_per_year(rows, "dra_cosdec_arcsec"),
        },
        "ddec_arcsec": {
            "median": median(de),
            "mad": mad(de),
            "rms": rms(de),
            "slope_arcsec_per_year": linear_slope_per_year(rows, "ddec_arcsec"),
        },
        "station_catalog_demeaned": {
            "dra_cosdec_rms_arcsec": rms(ra_dm),
            "ddec_rms_arcsec": rms(de_dm),
            "dra_cosdec_slope_arcsec_per_year": linear_slope_per_year(work, "dra_cosdec_demeaned_arcsec"),
            "ddec_slope_arcsec_per_year": linear_slope_per_year(work, "ddec_demeaned_arcsec"),
        },
    }


def object_run(spec):
    response = requests.get(
        MPC_URL,
        json={"desigs": [spec["designation"]], "output_format": ["ADES_DF"], "ades_version": "2022"},
        timeout=90,
    )
    response.raise_for_status()
    payload = response.json()
    raw_rows = payload[0]["ADES_DF"]

    usable = []
    skipped = collections.Counter()
    for row in raw_rows:
        ra = ffloat(row.get("ra"))
        dec = ffloat(row.get("dec"))
        obstime = row.get("obstime")
        stn = row.get("stn")
        if ra is None or dec is None:
            skipped["missing_ra_dec"] += 1
            continue
        if not obstime:
            skipped["missing_time"] += 1
            continue
        if not stn:
            skipped["missing_station"] += 1
            continue
        usable.append({
            "obsid": row.get("obsid"),
            "obstime": str(obstime),
            "stn": str(stn),
            "astcat": row.get("astcat"),
            "mode": row.get("mode"),
            "ra_obs_deg": ra,
            "dec_obs_deg": dec,
        })

    by_station = collections.defaultdict(list)
    for row in usable:
        by_station[row["stn"]].append(row)

    residuals = []
    failures = []
    station_meta = {}
    for stn in sorted(by_station):
        station_rows = sorted(by_station[stn], key=lambda r: r["obstime"])
        for part in chunks(station_rows, CHUNK):
            try:
                pred, meta = horizons_batch(spec["command"], stn, part)
                station_meta.setdefault(stn, meta)
                for obs, ep in zip(part, pred):
                    dra_deg = wrap_deg(obs["ra_obs_deg"] - ep["ra_deg"])
                    cosd = math.cos(math.radians(ep["dec_deg"]))
                    residuals.append({
                        **obs,
                        "jd": ep["jd"],
                        "ra_pred_deg": ep["ra_deg"],
                        "dec_pred_deg": ep["dec_deg"],
                        "dra_deg": dra_deg,
                        "dra_cosdec_arcsec": dra_deg * cosd * 3600.0,
                        "ddec_arcsec": (obs["dec_obs_deg"] - ep["dec_deg"]) * 3600.0,
                    })
            except Exception as exc:
                failures.append({
                    "station": stn,
                    "n_rows": len(part),
                    "first_obstime": part[0]["obstime"],
                    "last_obstime": part[-1]["obstime"],
                    "error_type": type(exc).__name__,
                    "error": str(exc),
                })
            time.sleep(0.08)

    residuals.sort(key=lambda r: r["jd"])
    canonical = json.dumps(residuals, sort_keys=True, separators=(",", ":")).encode("utf-8")
    years = collections.defaultdict(list)
    for r in residuals:
        year = int(r["obstime"][:4])
        years[year].append(r)

    return {
        "designation": spec["designation"],
        "command": spec["command"],
        "status": "PASS" if residuals and not failures else ("PARTIAL" if residuals else "FAIL"),
        "n_mpc_rows": len(raw_rows),
        "n_usable": len(usable),
        "n_residual_rows": len(residuals),
        "skipped": dict(skipped),
        "failures": failures,
        "first_obstime": residuals[0]["obstime"] if residuals else None,
        "last_obstime": residuals[-1]["obstime"] if residuals else None,
        "residual_sha256": hashlib.sha256(canonical).hexdigest(),
        "summary": summarize_residuals(residuals) if residuals else None,
        "yearly": {
            str(y): summarize_residuals(rs)
            for y, rs in sorted(years.items())
        },
        "top_stations": collections.Counter(r["stn"] for r in residuals).most_common(20),
        "top_catalogs": collections.Counter(r.get("astcat") or "" for r in residuals).most_common(20),
        "station_meta": station_meta,
        "largest_abs_residuals": sorted(
            residuals,
            key=lambda r: math.hypot(r["dra_cosdec_arcsec"], r["ddec_arcsec"]),
            reverse=True,
        )[:12],
    }


def main():
    results = []
    for spec in OBJECTS:
        r = object_run(spec)
        results.append(r)
        print("PBH_JPL_OC_OBJECT=" + json.dumps(r, sort_keys=True), flush=True)

    receipt = {
        "schema": "PBH_MPC_HORIZONS_OC_V0_2",
        "status": "PASS" if all(r["status"] == "PASS" for r in results) else "PARTIAL",
        "authority": "OBSERVED_MINUS_JPL_MODEL_ONLY",
        "mpc_endpoint": MPC_URL,
        "horizons_endpoint": HORIZONS_URL,
        "horizons_settings": {
            "EPHEM_TYPE": "OBSERVER",
            "QUANTITIES": "1",
            "ANG_FORMAT": "DEG",
            "CAL_FORMAT": "JD",
            "TIME_TYPE": "UT",
            "APPARENT": "AIRLESS",
            "force_model_as_reported_by_horizons": "DE441 + SB441-N16",
        },
        "results": results,
        "interpretation_boundary": {
            "jpl_orbit_is_fit_to_mpc_data": True,
            "residuals_are_not_independent_holdout_errors": True,
            "station_catalog_bias_not_yet_fully_modeled": True,
            "local_orbit_nuisance_not_yet_profiled": True,
            "common_tidal_tensor_not_yet_fitted": True,
            "PBH_not_inferred": True,
        },
        "canon_allowed": False,
    }
    print("PBH_JPL_OC_RECEIPT=" + json.dumps(receipt, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
