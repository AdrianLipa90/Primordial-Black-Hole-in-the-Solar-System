#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import math
import os
import re
import tempfile
from pathlib import Path

import requests
import empyrean

MPC_URL = "https://data.minorplanetcenter.net/api/get-obs"
OBJECTS = ["Sedna", "2012 VP113", "2015 TG387"]


def finite_or_none(x):
    try:
        v = float(x)
    except Exception:
        return None
    return v if math.isfinite(v) else None


def scalar(obj, name):
    try:
        return finite_or_none(getattr(obj, name))
    except Exception:
        return None


def sanitize(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", name)


def dataframe_numeric_summary(df):
    out = {}
    for col in df.columns:
        key = str(col)
        low = key.lower()
        if not any(k in low for k in ("resid", "chi", "weight", "cov", "mjd", "epoch")):
            continue
        try:
            vals = df[col].astype(float)
            vals = vals[vals.map(math.isfinite)]
        except Exception:
            continue
        if len(vals) == 0:
            continue
        out[key] = {
            "n": int(len(vals)),
            "mean": float(vals.mean()),
            "std": float(vals.std(ddof=1)) if len(vals) > 1 else 0.0,
            "min": float(vals.min()),
            "max": float(vals.max()),
        }
    return out


def fit_one(desig: str, root: Path) -> dict:
    resp = requests.get(
        MPC_URL,
        json={"desigs": [desig], "output_format": ["XML"], "ades_version": "2022"},
        timeout=90,
    )
    resp.raise_for_status()
    payload = resp.json()
    if not payload or "XML" not in payload[0]:
        raise RuntimeError(f"{desig}: MPC response lacks XML")
    xml = payload[0]["XML"]
    raw = xml.encode("utf-8")
    xml_path = root / f"{sanitize(desig)}.xml"
    xml_path.write_bytes(raw)

    obs, radar = empyrean.read_ades(str(xml_path))
    fits = empyrean.determine(obs, radar=radar)
    fit = fits.single()

    summary = fit.summary
    residuals = fit.observations
    rdf = residuals.to_dataframe()

    rejected = None
    try:
        rejected = int(len(residuals.rejected_only()))
    except Exception:
        pass

    return {
        "designation": desig,
        "mpc_xml_sha256": hashlib.sha256(raw).hexdigest(),
        "mpc_xml_bytes": len(raw),
        "n_input_optical": int(len(obs)),
        "n_input_radar": int(len(radar)),
        "converged": bool(fit.converged),
        "termination": str(fit.termination),
        "accepted_steps": getattr(fit, "accepted_steps", None),
        "final_solve_iterations": getattr(fit, "final_solve_iterations", None),
        "covariance_trust": str(getattr(fit, "covariance_trust", None)),
        "summary": {
            "reduced_chi2": scalar(summary, "reduced_chi2"),
            "rms_combined_arcsec": scalar(summary, "rms_combined_arcsec"),
            "rms_ra_arcsec": scalar(summary, "rms_ra_arcsec"),
            "rms_dec_arcsec": scalar(summary, "rms_dec_arcsec"),
        },
        "n_residual_rows": int(len(residuals)),
        "n_rejected": rejected,
        "residual_columns": [str(c) for c in rdf.columns],
        "numeric_residual_summary": dataframe_numeric_summary(rdf),
        "residual_head": json.loads(rdf.head(3).to_json(orient="records", date_format="iso")),
    }


def main() -> None:
    empyrean.initialize()
    root = Path(tempfile.mkdtemp(prefix="pbh-od-v02-"))
    results = []
    failures = []

    for desig in OBJECTS:
        try:
            result = fit_one(desig, root)
            results.append(result)
            print("PBH_OD_OBJECT=" + json.dumps(result, sort_keys=True), flush=True)
        except Exception as exc:
            failure = {
                "designation": desig,
                "error_type": type(exc).__name__,
                "error": str(exc),
            }
            failures.append(failure)
            print("PBH_OD_FAILURE=" + json.dumps(failure, sort_keys=True), flush=True)

    receipt = {
        "schema": "PBH_MPC_EMPYREAN_OD_PILOT_V0_2",
        "status": "PASS" if results and not failures else ("PARTIAL" if results else "FAIL"),
        "authority": "OBSERVATIONAL_RESIDUALS_ONLY",
        "objects_requested": OBJECTS,
        "n_delivered": len(results),
        "n_failed": len(failures),
        "results": results,
        "failures": failures,
        "mpc_endpoint": MPC_URL,
        "empyrean_version": empyrean.version_string(),
        "interpretation_boundary": {
            "standard_force_model_only": True,
            "external_tide_not_fitted_yet": True,
            "PBH_not_inferred": True,
            "purpose": "validate real MPC ADES -> independent N-body OD -> post-fit residual pipeline",
        },
        "canon_allowed": False,
    }
    print("PBH_OD_PILOT_RECEIPT=" + json.dumps(receipt, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
