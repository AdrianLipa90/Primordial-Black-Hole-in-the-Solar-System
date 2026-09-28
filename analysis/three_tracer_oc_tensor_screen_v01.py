#!/usr/bin/env python3
from __future__ import annotations

import csv
import io
import json
import math

import numpy as np
import requests

HORIZONS_URL = "https://ssd.jpl.nasa.gov/api/horizons.api"
ARCSEC_TO_RAD = math.pi / (180.0 * 3600.0)

# Frozen from PBH_MPC_HORIZONS_OC_RESIDUALS_V0_2
SOURCE_ROWS_SHA256 = "78492d3031d0e3fd41c883a14117c79c66e5c165d072b2c0293749c12d9ebbf3"

PROFILES = {
    "all_station_catalog_profiled": {
        "2012 VP113": {
            "command": "'DES=2012 VP113;'",
            "t0_jd_ut": 2457126.7521158545,
            "dra_arcsec_yr2": -1.1786631749570087e-05,
            "ddec_arcsec_yr2": 2.2615152683874284e-04,
        },
        "2015 TG387": {
            "command": "'541132;'",
            "t0_jd_ut": 2458704.6562090847,
            "dra_arcsec_yr2": -2.6247696089711193e-05,
            "ddec_arcsec_yr2": -1.2479664796914840e-05,
        },
        "Sedna": {
            "command": "'90377;'",
            "t0_jd_ut": 2458049.6644799025,
            "dra_arcsec_yr2": 6.806180461773833e-05,
            "ddec_arcsec_yr2": -1.9622805099622816e-04,
        },
    },
    "gaia_family_station_catalog_profiled": {
        "2012 VP113": {
            "command": "'DES=2012 VP113;'",
            "t0_jd_ut": 2457746.7402043464,
            "dra_arcsec_yr2": -6.961334931244616e-04,
            "ddec_arcsec_yr2": 1.8428030163705406e-04,
        },
        "2015 TG387": {
            "command": "'541132;'",
            "t0_jd_ut": 2458214.790899996,
            "dra_arcsec_yr2": -3.422434771802735e-05,
            "ddec_arcsec_yr2": -4.510116347653832e-05,
        },
        "Sedna": {
            "command": "'90377;'",
            "t0_jd_ut": 2459070.5574593414,
            "dra_arcsec_yr2": -1.9724662276312152e-04,
            "ddec_arcsec_yr2": -7.370083810209547e-04,
        },
    },
}


def horizons_radec(command, jd_ut):
    params = {
        "format": "json",
        "COMMAND": command,
        "OBJ_DATA": "'NO'",
        "MAKE_EPHEM": "'YES'",
        "EPHEM_TYPE": "'OBSERVER'",
        "CENTER": "'500@399'",
        "TLIST": f"'{jd_ut:.12f}'",
        "TLIST_TYPE": "'JD'",
        "TIME_TYPE": "'UT'",
        "QUANTITIES": "'1'",
        "REF_SYSTEM": "'ICRF'",
        "ANG_FORMAT": "'DEG'",
        "CSV_FORMAT": "'YES'",
        "EXTRA_PREC": "'YES'",
        "CAL_FORMAT": "'BOTH'",
    }
    r = requests.get(HORIZONS_URL, params=params, timeout=60)
    r.raise_for_status()
    result = r.json().get("result", "")
    if "** Error" in result or "$$SOE" not in result:
        raise RuntimeError(result[:1000])
    body = result.split("$$SOE", 1)[1].split("$$EOE", 1)[0].strip()
    rec = next(csv.reader(io.StringIO(body)))
    vals = [x.strip() for x in rec]
    return float(vals[4]), float(vals[5])


def basis(ra_deg, dec_deg):
    a = math.radians(ra_deg)
    d = math.radians(dec_deg)
    n = np.array([math.cos(d)*math.cos(a), math.cos(d)*math.sin(a), math.sin(d)])
    era = np.array([-math.sin(a), math.cos(a), 0.0])
    edec = np.array([-math.cos(a)*math.sin(d), -math.sin(a)*math.sin(d), math.cos(d)])
    return n, era, edec


def design_row(e, n):
    # T = [[t1,t2,t3],[t2,t4,t5],[t3,t5,-t1-t4]]
    return np.array([
        e[0]*n[0] - e[2]*n[2],
        e[0]*n[1] + e[1]*n[0],
        e[0]*n[2] + e[2]*n[0],
        e[1]*n[1] - e[2]*n[2],
        e[1]*n[2] + e[2]*n[1],
    ])


def vec5_to_tensor(t):
    return np.array([
        [t[0], t[1], t[2]],
        [t[1], t[3], t[4]],
        [t[2], t[4], -t[0]-t[3]],
    ])


def point_mass_diagnostic(T):
    vals, vecs = np.linalg.eigh(T)
    axis = vecs[:, -1]
    # deterministic sign only; physical sign remains degenerate.
    for x in axis:
        if abs(x) > 1e-15:
            if x < 0:
                axis = -axis
            break
    Q = 3.0*np.outer(axis, axis) - np.eye(3)
    mu = float(np.sum(T*Q)/6.0)
    P = mu*Q
    shape = float(np.linalg.norm(T-P, "fro")/max(np.linalg.norm(T, "fro"), 1e-300))
    split = float(abs(vals[1]-vals[0])/max(abs(mu), 1e-300))
    passes = {
        "positive_mu": mu > 0,
        "point_mass_shape": shape <= 0.20,
        "transverse_degeneracy": split <= 0.40,
    }
    passes["all"] = all(passes.values())
    return {
        "eigenvalues_per_year2": vals.tolist(),
        "axis_sign_degenerate": axis.tolist(),
        "mu_per_year2": mu,
        "point_mass_shape_residual": shape,
        "transverse_eigenvalue_split_over_mu": split,
        "thresholds": {"max_shape_residual": 0.20, "max_transverse_split_over_mu": 0.40},
        "passes": passes,
    }


def fit_profile(name, rows):
    H = []
    y = []
    geometry = {}
    for obj, d in rows.items():
        ra, dec = horizons_radec(d["command"], d["t0_jd_ut"])
        n, era, edec = basis(ra, dec)
        geometry[obj] = {
            "t0_jd_ut": d["t0_jd_ut"],
            "ra_deg": ra,
            "dec_deg": dec,
            "n_icrf": n.tolist(),
            "curvature_arcsec_yr2": [d["dra_arcsec_yr2"], d["ddec_arcsec_yr2"]],
        }
        H.append(design_row(era, n))
        H.append(design_row(edec, n))
        y.append(d["dra_arcsec_yr2"] * ARCSEC_TO_RAD)
        y.append(d["ddec_arcsec_yr2"] * ARCSEC_TO_RAD)

    H = np.asarray(H, dtype=float)
    y = np.asarray(y, dtype=float)
    t, residuals, rank, singular = np.linalg.lstsq(H, y, rcond=None)
    T = vec5_to_tensor(t)
    pred = H @ t
    resid = y - pred
    ynorm = float(np.linalg.norm(y))
    return {
        "name": name,
        "geometry": geometry,
        "rank": int(rank),
        "singular_values": singular.tolist(),
        "condition_number": float(singular[0]/singular[-1]) if singular[-1] > 0 else float("inf"),
        "tensor_vec5_per_year2": t.tolist(),
        "tensor_3x3_per_year2": T.tolist(),
        "measurement_norm_rad_per_year2": ynorm,
        "fit_residual_norm_rad_per_year2": float(np.linalg.norm(resid)),
        "relative_fit_residual": float(np.linalg.norm(resid)/ynorm) if ynorm > 0 else None,
        "predicted_measurements_rad_per_year2": pred.tolist(),
        "observed_measurements_rad_per_year2": y.tolist(),
        "point_mass": point_mass_diagnostic(T),
    }


def frob_cos(A, B):
    a=np.asarray(A); b=np.asarray(B)
    return float(np.sum(a*b)/(np.linalg.norm(a,"fro")*np.linalg.norm(b,"fro")))


def main():
    fits = {name: fit_profile(name, p) for name, p in PROFILES.items()}
    a = fits["all_station_catalog_profiled"]
    g = fits["gaia_family_station_catalog_profiled"]
    axis_overlap = abs(float(np.dot(
        a["point_mass"]["axis_sign_degenerate"],
        g["point_mass"]["axis_sign_degenerate"],
    )))
    comparison = {
        "axis_overlap_abs_all_vs_gaia": axis_overlap,
        "tensor_frobenius_cosine_all_vs_gaia": frob_cos(
            a["tensor_3x3_per_year2"], g["tensor_3x3_per_year2"]
        ),
        "both_point_mass_shape_pass": bool(
            a["point_mass"]["passes"]["all"] and g["point_mass"]["passes"]["all"]
        ),
        "stability_gate": bool(
            axis_overlap >= 0.95
            and abs(frob_cos(a["tensor_3x3_per_year2"], g["tensor_3x3_per_year2"])) >= 0.90
        ),
    }
    receipt = {
        "schema": "PBH_THREE_TRACER_OC_TENSOR_SCREEN_V0_1",
        "status": "SCREEN_ONLY",
        "source_rows_sha256": SOURCE_ROWS_SHA256,
        "fits": fits,
        "comparison": comparison,
        "interpretation_boundary": {
            "input_curvatures_are_postfit_OC_proxies": True,
            "current_JPL_solution_not_independent_of_MPC_astrometry": True,
            "do_not_convert_mu_to_mass_distance": True,
            "not_a_variational_orbit_fit": True,
            "PBH_not_inferred": True,
        },
        "canon_allowed": False,
    }
    print("PBH_OC_TENSOR_SCREEN=" + json.dumps(receipt, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
