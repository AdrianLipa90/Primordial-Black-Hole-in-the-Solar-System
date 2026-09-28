#!/usr/bin/env python3
from __future__ import annotations

import csv
import hashlib
import io
import json
import math
import statistics
from collections import Counter, defaultdict
from datetime import datetime, timezone

import requests
import numpy as np

MPC_URL = "https://data.minorplanetcenter.net/api/get-obs"
HORIZONS_URL = "https://ssd.jpl.nasa.gov/api/horizons.api"
CHUNK = 40

OBJECTS = {
    "2012 VP113": "'DES=2012 VP113;'",
    "2015 TG387": "'541132;'",
    "Sedna": "'90377;'",
}


def as_float(x):
    try:
        v = float(x)
        return v if math.isfinite(v) else None
    except Exception:
        return None


def parse_iso_utc(s):
    return datetime.fromisoformat(s.replace("Z", "+00:00")).astimezone(timezone.utc)


def horizons_time(s):
    dt = parse_iso_utc(s)
    return dt.strftime("%Y-%b-%d %H:%M:%S.") + f"{dt.microsecond // 1000:03d}"


def fetch_mpc(desig):
    r = requests.get(
        MPC_URL,
        json={"desigs": [desig], "output_format": ["ADES_DF"], "ades_version": "2022"},
        timeout=90,
    )
    r.raise_for_status()
    payload = r.json()
    if not payload or "ADES_DF" not in payload[0]:
        raise RuntimeError(f"{desig}: MPC ADES_DF missing")
    out = []
    for row in payload[0]["ADES_DF"]:
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
        })
    return out


def parse_horizons_result(result):
    if "$$SOE" not in result or "$$EOE" not in result:
        raise RuntimeError("Horizons result has no SOE/EOE block")
    body = result.split("$$SOE", 1)[1].split("$$EOE", 1)[0].strip()
    rows = []
    for rec in csv.reader(io.StringIO(body)):
        vals = [x.strip() for x in rec]
        if len(vals) < 6:
            continue
        jd = as_float(vals[1])
        ra = as_float(vals[4])
        dec = as_float(vals[5])
        if jd is None or ra is None or dec is None:
            raise RuntimeError(f"unparseable Horizons row: {vals!r}")
        rows.append({
            "cal": vals[0],
            "jd_ut": jd,
            "ra_pred_deg": ra,
            "dec_pred_deg": dec,
        })
    return rows


def query_horizons_chunk(command, stn, obs_rows):
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
        raise RuntimeError(f"Horizons error for {stn}: {result[:700]}")
    pred = parse_horizons_result(result)
    if len(pred) != len(obs_rows):
        raise RuntimeError(
            f"{stn}: Horizons returned {len(pred)} rows for {len(obs_rows)} observations"
        )
    return pred


def query_horizons(command, stn, obs_rows):
    out = []
    for i in range(0, len(obs_rows), CHUNK):
        out.extend(query_horizons_chunk(command, stn, obs_rows[i:i+CHUNK]))
    if len(out) != len(obs_rows):
        raise RuntimeError(f"{stn}: chunk assembly length mismatch")
    return out


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
        "ra_obs_deg": obs["ra_obs_deg"],
        "dec_obs_deg": obs["dec_obs_deg"],
        "ra_pred_deg": pred["ra_pred_deg"],
        "dec_pred_deg": pred["dec_pred_deg"],
        "dra_cosdec_arcsec": dra_as,
        "ddec_arcsec": ddec_as,
        "rho_arcsec": math.hypot(dra_as, ddec_as),
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


def summary(rows):
    return {
        "dra": med_mad([r["dra_cosdec_arcsec"] for r in rows]),
        "ddec": med_mad([r["ddec_arcsec"] for r in rows]),
        "rho": med_mad([r["rho_arcsec"] for r in rows]),
    }


def solve3(a, b):
    m = [list(map(float, row)) + [float(b[i])] for i, row in enumerate(a)]
    for col in range(3):
        piv = max(range(col, 3), key=lambda i: abs(m[i][col]))
        if abs(m[piv][col]) < 1e-20:
            return None
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


def quadratic_fit_values(rows, values):
    if len(rows) < 5:
        return None
    t0 = sum(r["jd_ut"] for r in rows) / len(rows)
    ts = [(r["jd_ut"] - t0) / 365.25 for r in rows]
    ys = list(map(float, values))
    s0 = len(ts)
    s1 = sum(ts)
    s2 = sum(t*t for t in ts)
    s3 = sum(t**3 for t in ts)
    s4 = sum(t**4 for t in ts)
    b0 = sum(ys)
    b1 = sum(t*y for t, y in zip(ts, ys))
    b2 = sum(t*t*y for t, y in zip(ts, ys))
    c = solve3(
        [[s0, s1, s2], [s1, s2, s3], [s2, s3, s4]],
        [b0, b1, b2],
    )
    if c is None:
        return None
    resid = [y - (c[0] + c[1]*t + c[2]*t*t) for t, y in zip(ts, ys)]
    return {
        "n": len(rows),
        "t0_jd_ut": t0,
        "intercept_arcsec": c[0],
        "linear_arcsec_per_year": c[1],
        "quadratic_arcsec_per_year2": c[2],
        "second_derivative_arcsec_per_year2": 2.0*c[2],
        "postfit_rms_arcsec": math.sqrt(sum(x*x for x in resid) / len(resid)),
    }


def _nuisance_matrix(rows):
    """Nuisance design Z=[1,t,station dummies,catalog dummies].

    The quadratic regressor t^2 is intentionally excluded from Z.  FWL then
    residualizes *both* y and t^2 against exactly the same nuisance space.
    """
    n=len(rows)
    t0=sum(r["jd_ut"] for r in rows)/n
    t=np.asarray([(r["jd_ut"]-t0)/365.25 for r in rows],dtype=float)
    cols=[np.ones(n,dtype=float),t]
    labels=["intercept","linear_time"]

    for key in ("stn","catalog"):
        levels=sorted({str(r.get(key) or "UNKNOWN") for r in rows})
        # Drop one reference level. Any remaining rank deficiency is handled
        # by SVD pseudoinverse and reported explicitly.
        for level in levels[1:]:
            cols.append(np.asarray(
                [1.0 if str(r.get(key) or "UNKNOWN")==level else 0.0 for r in rows],
                dtype=float,
            ))
            labels.append(f"{key}:{level}")
    Z=np.column_stack(cols)
    return t0,t,Z,labels


def fwl_quadratic_fit(rows, field):
    if len(rows)<5:
        return None
    t0,t,Z,labels=_nuisance_matrix(rows)
    y=np.asarray([float(r[field]) for r in rows],dtype=float)
    q=t*t

    pinv=np.linalg.pinv(Z,rcond=1e-12)
    y_perp=y-Z@(pinv@y)
    q_perp=q-Z@(pinv@q)

    qq=float(q_perp@q_perp)
    if qq<=1e-24:
        return {
            "status":"UNIDENTIFIABLE",
            "n":len(rows),"t0_jd_ut":t0,
            "rank_nuisance":int(np.linalg.matrix_rank(Z)),
            "n_nuisance_columns":int(Z.shape[1]),
            "q_perp_norm2":qq,
        }

    c2=float(q_perp@y_perp/qq)
    resid=y_perp-c2*q_perp
    rank_z=int(np.linalg.matrix_rank(Z))
    dof=max(1,len(rows)-rank_z-1)
    sigma2=float(resid@resid/dof)
    se_c2=math.sqrt(max(0.0,sigma2/qq))
    centered_q=q-q.mean()
    ident_frac=qq/max(float(centered_q@centered_q),1e-300)

    return {
        "status":"IDENTIFIED",
        "n":len(rows),
        "t0_jd_ut":t0,
        "quadratic_arcsec_per_year2":c2,
        "second_derivative_arcsec_per_year2":2.0*c2,
        "naive_ols_se_second_arcsec_per_year2":2.0*se_c2,
        "postfit_rms_arcsec":math.sqrt(float(resid@resid)/len(rows)),
        "rank_nuisance":rank_z,
        "n_nuisance_columns":int(Z.shape[1]),
        "dof_resid":dof,
        "q_perp_norm2":qq,
        "quadratic_identifiability_fraction_after_nuisance":ident_frac,
        "note":"SE is iid-OLS diagnostic only; not a calibrated significance.",
    }


def curvature(rows, profiled=False):
    if len(rows)<5:
        return None
    out={}
    comps=[]
    for field,key in (("dra_cosdec_arcsec","dra"),("ddec_arcsec","ddec")):
        fit=fwl_quadratic_fit(rows,field) if profiled else quadratic_fit_values(
            rows,[r[field] for r in rows]
        )
        out[key]=fit
        if fit is None or fit.get("status")=="UNIDENTIFIABLE":
            return None
        comps.append(fit["second_derivative_arcsec_per_year2"])
    out["vector_norm_arcsec_per_year2"]=math.hypot(*comps)
    out["profiling"]="FWL_joint_station_catalog" if profiled else "raw_quadratic"
    return out


def range_summary(values):
    vals = sorted(v for v in values if v is not None and math.isfinite(v))
    if not vals:
        return None
    return {
        "n": len(vals),
        "min": vals[0],
        "median": vals[len(vals)//2],
        "max": vals[-1],
    }


def leave_one_group(rows, key):
    levels = sorted({str(r.get(key) or "UNKNOWN") for r in rows})
    vectors = []
    for level in levels:
        sub = [r for r in rows if str(r.get(key) or "UNKNOWN") != level]
        if len(sub) < 8:
            continue
        c = curvature(sub, profiled=True)
        if c is not None:
            vectors.append({
                "removed": level,
                "n": len(sub),
                "norm": c["vector_norm_arcsec_per_year2"],
                "dra": c["dra"]["second_derivative_arcsec_per_year2"],
                "ddec": c["ddec"]["second_derivative_arcsec_per_year2"],
            })
    return {
        "n_runs": len(vectors),
        "norm_range": range_summary([x["norm"] for x in vectors]),
        "dra_range": range_summary([x["dra"] for x in vectors]),
        "ddec_range": range_summary([x["ddec"] for x in vectors]),
        "largest_norm_cases": sorted(vectors, key=lambda x: x["norm"], reverse=True)[:3],
        "runs": vectors,
    }


def top_group_offsets(rows, key, top=8):
    groups = defaultdict(list)
    for r in rows:
        groups[str(r.get(key) or "UNKNOWN")].append(r)
    chosen = sorted(groups.items(), key=lambda kv: len(kv[1]), reverse=True)[:top]
    return {
        k: {
            "n": len(rr),
            "mean_dra_arcsec": sum(r["dra_cosdec_arcsec"] for r in rr)/len(rr),
            "mean_ddec_arcsec": sum(r["ddec_arcsec"] for r in rr)/len(rr),
            "rho_rms_arcsec": math.sqrt(sum(r["rho_arcsec"]**2 for r in rr)/len(rr)),
        }
        for k, rr in chosen
    }


def object_diagnostic(rows, n_mpc_usable, station_failures):
    raw=curvature(rows,profiled=False)
    prof=curvature(rows,profiled=True)
    gaia=[r for r in rows if str(r.get("catalog") or "").startswith("Gaia")]
    gprof=curvature(gaia,profiled=True) if len(gaia)>=5 else None
    suppression=(
        raw["vector_norm_arcsec_per_year2"]/prof["vector_norm_arcsec_per_year2"]
        if raw and prof and prof["vector_norm_arcsec_per_year2"]>0 else None
    )
    return {
        "n_mpc_usable":n_mpc_usable,
        "n_residuals":len(rows),
        "n_station_failures":len(station_failures),
        "station_failures":station_failures,
        "raw_summary":summary(rows),
        "raw_curvature":raw,
        "station_catalog_profiled_curvature":prof,
        "curvature_norm_suppression_factor":suppression,
        "gaia_family_n":len(gaia),
        "gaia_family_profiled_curvature":gprof,
        "top_station_offsets":top_group_offsets(rows,"stn"),
        "top_catalog_offsets":top_group_offsets(rows,"catalog"),
        "stability_layer":"DEFERRED_TO_OPTIMIZED_LOSO",
    }

def main():
    all_rows = []
    objects = {}
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
                print("PBH_OC_V04_PROGRESS="+json.dumps({"object":obj,"station":stn,"n":len(rr)},sort_keys=True),flush=True)
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
            objects[obj] = object_diagnostic(object_rows, len(obs), station_failures)
        except Exception as exc:
            failures.append({
                "object": obj,
                "error_type": type(exc).__name__,
                "error": str(exc),
            })

    canonical = json.dumps(all_rows, sort_keys=True, separators=(",", ":")).encode()
    receipt = {
        "schema": "PBH_MPC_HORIZONS_OC_RESIDUALS_V0_4_CORE_FWL",
        "script_revision": "three-tracer-full-jackknife-r2",
        "status": (
            "PASS" if len(objects) == len(OBJECTS) and not failures
            and all(v["n_station_failures"] == 0 for v in objects.values())
            else ("PARTIAL" if objects else "FAIL")
        ),
        "authority": "OBSERVATIONAL_OC_DIAGNOSTIC_ONLY",
        "mpc_endpoint": MPC_URL,
        "horizons_endpoint": HORIZONS_URL,
        "horizons_baseline": "current JPL solution; not statistically independent of MPC astrometry",
        "objects": objects,
        "n_rows": len(all_rows),
        "rows_sha256": hashlib.sha256(canonical).hexdigest(),
        "failures": failures,
        "interpretation_boundary": {
            "not_a_PBH_detection": True,
            "not_independent_of_JPL_orbit_fit": True,
            "polynomial_curvature_not_a_tidal_tensor": True,
            "purpose": "three-tracer O-C systematics screen with correct FWL nuisance profiling before variational common-tide inversion",
        },
        "canon_allowed": False,
    }
    compact={
        "schema":"PBH_MPC_HORIZONS_OC_FWL_SUMMARY_V0_4",
        "status":receipt["status"],
        "rows_sha256":receipt["rows_sha256"],
        "n_rows":receipt["n_rows"],
        "objects":{},
    }
    for name,o in objects.items():
        compact["objects"][name]={
            "n_residuals":o["n_residuals"],
            "n_station_failures":o["n_station_failures"],
            "raw_vector_norm":None if o["raw_curvature"] is None else o["raw_curvature"]["vector_norm_arcsec_per_year2"],
            "fwl_vector_norm":None if o["station_catalog_profiled_curvature"] is None else o["station_catalog_profiled_curvature"]["vector_norm_arcsec_per_year2"],
            "fwl_dra_second":None if o["station_catalog_profiled_curvature"] is None else o["station_catalog_profiled_curvature"]["dra"]["second_derivative_arcsec_per_year2"],
            "fwl_ddec_second":None if o["station_catalog_profiled_curvature"] is None else o["station_catalog_profiled_curvature"]["ddec"]["second_derivative_arcsec_per_year2"],
            "fwl_dra_se_naive":None if o["station_catalog_profiled_curvature"] is None else o["station_catalog_profiled_curvature"]["dra"]["naive_ols_se_second_arcsec_per_year2"],
            "fwl_ddec_se_naive":None if o["station_catalog_profiled_curvature"] is None else o["station_catalog_profiled_curvature"]["ddec"]["naive_ols_se_second_arcsec_per_year2"],
            "quadratic_identifiability_dra":None if o["station_catalog_profiled_curvature"] is None else o["station_catalog_profiled_curvature"]["dra"]["quadratic_identifiability_fraction_after_nuisance"],
            "quadratic_identifiability_ddec":None if o["station_catalog_profiled_curvature"] is None else o["station_catalog_profiled_curvature"]["ddec"]["quadratic_identifiability_fraction_after_nuisance"],
            "leave_one_station_norm_range":o["leave_one_station"]["norm_range"],
            "leave_one_catalog_norm_range":o["leave_one_catalog"]["norm_range"],
        }
    print("PBH_OC_V04_SUMMARY="+json.dumps(compact,sort_keys=True),flush=True)
    print("PBH_OC_V04_RECEIPT=" + json.dumps(receipt, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
