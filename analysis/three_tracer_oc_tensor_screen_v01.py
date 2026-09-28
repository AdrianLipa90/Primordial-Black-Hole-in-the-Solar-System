#!/usr/bin/env python3
from __future__ import annotations

import csv
import io
import json
import math

import requests

HORIZONS_URL = "https://ssd.jpl.nasa.gov/api/horizons.api"
ARCSEC_TO_RAD = math.pi / (180.0 * 3600.0)
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
    vals = [x.strip() for x in next(csv.reader(io.StringIO(body)))]
    return float(vals[4]), float(vals[5])


def dot(a, b):
    return sum(x*y for x, y in zip(a, b))


def norm(v):
    return math.sqrt(dot(v, v))


def matvec(A, x):
    return [dot(row, x) for row in A]


def basis(ra_deg, dec_deg):
    a = math.radians(ra_deg)
    d = math.radians(dec_deg)
    n = [math.cos(d)*math.cos(a), math.cos(d)*math.sin(a), math.sin(d)]
    era = [-math.sin(a), math.cos(a), 0.0]
    edec = [-math.cos(a)*math.sin(d), -math.sin(a)*math.sin(d), math.cos(d)]
    return n, era, edec


def design_row(e, n):
    return [
        e[0]*n[0] - e[2]*n[2],
        e[0]*n[1] + e[1]*n[0],
        e[0]*n[2] + e[2]*n[0],
        e[1]*n[1] - e[2]*n[2],
        e[1]*n[2] + e[2]*n[1],
    ]


def solve_linear(A, b):
    n = len(b)
    m = [list(map(float, A[i])) + [float(b[i])] for i in range(n)]
    pivots = []
    for col in range(n):
        p = max(range(col, n), key=lambda i: abs(m[i][col]))
        if abs(m[p][col]) < 1e-30:
            raise RuntimeError("singular normal matrix")
        m[col], m[p] = m[p], m[col]
        pivots.append(abs(m[col][col]))
        q = m[col][col]
        for j in range(col, n+1):
            m[col][j] /= q
        for i in range(n):
            if i == col:
                continue
            q = m[i][col]
            for j in range(col, n+1):
                m[i][j] -= q*m[col][j]
    return [m[i][n] for i in range(n)], pivots


def normal_equations(H, y):
    p = len(H[0])
    A = [[0.0]*p for _ in range(p)]
    b = [0.0]*p
    for row, yi in zip(H, y):
        for i in range(p):
            b[i] += row[i]*yi
            for j in range(p):
                A[i][j] += row[i]*row[j]
    return A, b


def jacobi_eigh(A, iterations=256, tol=1e-16):
    n = len(A)
    a = [list(map(float, r)) for r in A]
    V = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    for _ in range(iterations):
        p, q, largest = 0, 1, 0.0
        for i in range(n):
            for j in range(i+1, n):
                z = abs(a[i][j])
                if z > largest:
                    p, q, largest = i, j, z
        scale = max(1.0, max(abs(x) for row in a for x in row))
        if largest <= tol*scale:
            break
        phi = 0.5*math.atan2(2.0*a[p][q], a[q][q]-a[p][p])
        c, s = math.cos(phi), math.sin(phi)
        app, aqq, apq = a[p][p], a[q][q], a[p][q]
        for k in range(n):
            if k in (p, q):
                continue
            akp, akq = a[k][p], a[k][q]
            a[k][p] = a[p][k] = c*akp - s*akq
            a[k][q] = a[q][k] = s*akp + c*akq
        a[p][p] = c*c*app - 2*s*c*apq + s*s*aqq
        a[q][q] = s*s*app + 2*s*c*apq + c*c*aqq
        a[p][q] = a[q][p] = 0.0
        for k in range(n):
            vkp, vkq = V[k][p], V[k][q]
            V[k][p] = c*vkp - s*vkq
            V[k][q] = s*vkp + c*vkq
    vals = [a[i][i] for i in range(n)]
    cols = [[V[r][i] for r in range(n)] for i in range(n)]
    order = sorted(range(n), key=lambda i: vals[i])
    return [vals[i] for i in order], [cols[i] for i in order]


def vec5_to_tensor(t):
    return [
        [t[0], t[1], t[2]],
        [t[1], t[3], t[4]],
        [t[2], t[4], -t[0]-t[3]],
    ]


def frob(A, B):
    return sum(A[i][j]*B[i][j] for i in range(3) for j in range(3))


def frob_norm(A):
    return math.sqrt(frob(A, A))


def point_mass_diagnostic(T):
    vals, vecs = jacobi_eigh(T, iterations=128)
    axis = list(vecs[-1])
    na = norm(axis)
    axis = [x/na for x in axis]
    for x in axis:
        if abs(x) > 1e-15:
            if x < 0:
                axis = [-y for y in axis]
            break
    Q = [[3.0*axis[i]*axis[j] - (1.0 if i == j else 0.0) for j in range(3)] for i in range(3)]
    mu = frob(T, Q)/6.0
    P = [[mu*Q[i][j] for j in range(3)] for i in range(3)]
    R = [[T[i][j]-P[i][j] for j in range(3)] for i in range(3)]
    shape = frob_norm(R)/max(frob_norm(T), 1e-300)
    split = abs(vals[1]-vals[0])/max(abs(mu), 1e-300)
    passes = {
        "positive_mu": mu > 0.0,
        "point_mass_shape": shape <= 0.20,
        "transverse_degeneracy": split <= 0.40,
    }
    passes["all"] = all(passes.values())
    return {
        "eigenvalues_per_year2": vals,
        "axis_sign_degenerate": axis,
        "mu_per_year2": mu,
        "point_mass_shape_residual": shape,
        "transverse_eigenvalue_split_over_mu": split,
        "thresholds": {"max_shape_residual": 0.20, "max_transverse_split_over_mu": 0.40},
        "passes": passes,
    }


def fit_profile(name, rows):
    H, y, geometry = [], [], {}
    for obj, d in rows.items():
        ra, dec = horizons_radec(d["command"], d["t0_jd_ut"])
        n, era, edec = basis(ra, dec)
        geometry[obj] = {
            "t0_jd_ut": d["t0_jd_ut"],
            "ra_deg": ra,
            "dec_deg": dec,
            "n_icrf": n,
            "curvature_arcsec_yr2": [d["dra_arcsec_yr2"], d["ddec_arcsec_yr2"]],
        }
        H.append(design_row(era, n))
        H.append(design_row(edec, n))
        y.extend([
            d["dra_arcsec_yr2"]*ARCSEC_TO_RAD,
            d["ddec_arcsec_yr2"]*ARCSEC_TO_RAD,
        ])

    N, g = normal_equations(H, y)
    t, pivots = solve_linear(N, g)
    T = vec5_to_tensor(t)
    pred = matvec(H, t)
    resid = [a-b for a, b in zip(y, pred)]
    y_norm = norm(y)

    gram_vals, _ = jacobi_eigh(N, iterations=512)
    positive = [max(v, 0.0) for v in gram_vals]
    svals = [math.sqrt(v) for v in positive]
    s_pos = [x for x in svals if x > 1e-15]
    condition = max(s_pos)/min(s_pos) if s_pos else float("inf")
    rank = len(s_pos)

    return {
        "name": name,
        "geometry": geometry,
        "rank": rank,
        "singular_values_approx": svals,
        "condition_number_approx": condition,
        "gaussian_pivot_ratio": max(pivots)/min(pivots),
        "tensor_vec5_per_year2": t,
        "tensor_3x3_per_year2": T,
        "measurement_norm_rad_per_year2": y_norm,
        "fit_residual_norm_rad_per_year2": norm(resid),
        "relative_fit_residual": norm(resid)/y_norm if y_norm > 0 else None,
        "predicted_measurements_rad_per_year2": pred,
        "observed_measurements_rad_per_year2": y,
        "point_mass": point_mass_diagnostic(T),
    }


def frob_cos(A, B):
    return frob(A, B)/(frob_norm(A)*frob_norm(B))


def main():
    fits = {name: fit_profile(name, p) for name, p in PROFILES.items()}
    a = fits["all_station_catalog_profiled"]
    g = fits["gaia_family_station_catalog_profiled"]
    axis_overlap = abs(dot(
        a["point_mass"]["axis_sign_degenerate"],
        g["point_mass"]["axis_sign_degenerate"],
    ))
    tensor_cos = frob_cos(a["tensor_3x3_per_year2"], g["tensor_3x3_per_year2"])
    comparison = {
        "axis_overlap_abs_all_vs_gaia": axis_overlap,
        "tensor_frobenius_cosine_all_vs_gaia": tensor_cos,
        "both_point_mass_shape_pass": bool(
            a["point_mass"]["passes"]["all"] and g["point_mass"]["passes"]["all"]
        ),
        "stability_gate": bool(axis_overlap >= 0.95 and abs(tensor_cos) >= 0.90),
    }
    receipt = {
        "schema": "PBH_THREE_TRACER_OC_TENSOR_SCREEN_V0_1",
        "script_revision": "pure-python-gauss-jacobi-r2",
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
