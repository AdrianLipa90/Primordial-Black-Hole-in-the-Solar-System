#!/usr/bin/env python3
from __future__ import annotations

import csv
import hashlib
import io
import json
import math
from collections import defaultdict
from datetime import datetime, timezone

import numpy as np
import rebound
import requests

MPC_URL = "https://data.minorplanetcenter.net/api/get-obs"
HORIZONS_URL = "https://ssd.jpl.nasa.gov/api/horizons.api"
EXPECTED_RAW_SHA = "78492d3031d0e3fd41c883a14117c79c66e5c165d072b2c0293749c12d9ebbf3"
DAY_PER_YEAR = 365.25
ARCSEC_PER_RAD = 180.0 * 3600.0 / math.pi
CHUNK = 40

DPOS_AU = 1.0e-5
DVEL_AU_YR = 1.0e-6
DTIDE_YR2 = 1.0e-9
RCOND = 1.0e-10

OBJECTS = {
    "2012 VP113": "'DES=2012 VP113;'",
    "2015 TG387": "'541132;'",
    "Sedna": "'90377;'",
}

# Approximate mass ratios to M_sun, adequate only for this derivative pilot.
# The states themselves are current Horizons states.
BODIES = [
    ("Mercury", "'199'", 1.660120e-7),
    ("Venus",   "'299'", 2.4478383e-6),
    ("Earth",   "'399'", 3.0034896e-6),
    ("Moon",    "'301'", 3.694303e-8),
    ("Mars",    "'499'", 3.2271514e-7),
    ("Jupiter", "'599'", 9.5479194e-4),
    ("Saturn",  "'699'", 2.8588598e-4),
    ("Uranus",  "'799'", 4.3662440e-5),
    ("Neptune", "'899'", 5.1513890e-5),
]

QBASIS = np.array([
    [[1.,0.,0.],[0.,0.,0.],[0.,0.,-1.]],
    [[0.,1.,0.],[1.,0.,0.],[0.,0.,0.]],
    [[0.,0.,1.],[0.,0.,0.],[1.,0.,0.]],
    [[0.,0.,0.],[0.,1.,0.],[0.,0.,-1.]],
    [[0.,0.,0.],[0.,0.,1.],[0.,1.,0.]],
])


def finite(x):
    try:
        v = float(x)
        return v if math.isfinite(v) else None
    except Exception:
        return None


def htime(s):
    dt = datetime.fromisoformat(s.replace("Z", "+00:00")).astimezone(timezone.utc)
    return dt.strftime("%Y-%b-%d %H:%M:%S.") + f"{dt.microsecond // 1000:03d}"


def fetch_mpc(desig):
    r = requests.get(
        MPC_URL,
        json={"desigs": [desig], "output_format": ["ADES_DF"], "ades_version": "2022"},
        timeout=90,
    )
    r.raise_for_status()
    p = r.json()
    if not p or "ADES_DF" not in p[0]:
        raise RuntimeError(f"{desig}: ADES_DF missing")
    out = []
    for row in p[0]["ADES_DF"]:
        ra, dec = finite(row.get("ra")), finite(row.get("dec"))
        stn = str(row.get("stn") or "").strip()
        tm = str(row.get("obstime") or "").strip()
        if ra is None or dec is None or not stn or not tm:
            continue
        out.append({
            "obs_id": row.get("obsid"),
            "time": tm,
            "stn": stn,
            "catalog": row.get("astcat"),
            "mode": row.get("mode"),
            "ra_obs_deg": ra,
            "dec_obs_deg": dec,
        })
    return out


def parse_observer(result):
    if "$$SOE" not in result or "$$EOE" not in result:
        raise RuntimeError("observer response missing SOE/EOE")
    body = result.split("$$SOE", 1)[1].split("$$EOE", 1)[0].strip()
    out = []
    for rec in csv.reader(io.StringIO(body)):
        v = [x.strip() for x in rec]
        if len(v) < 6:
            continue
        jd, ra, dec = finite(v[1]), finite(v[4]), finite(v[5])
        if jd is None or ra is None or dec is None:
            raise RuntimeError(f"bad observer row: {v!r}")
        out.append((jd, ra, dec))
    return out


def observer_query(command, stn, rows):
    out = []
    for i in range(0, len(rows), CHUNK):
        rr = rows[i:i+CHUNK]
        tlist = " ".join(f"'{htime(x['time'])}'" for x in rr)
        params = {
            "format": "json", "COMMAND": command, "OBJ_DATA": "'NO'",
            "MAKE_EPHEM": "'YES'", "EPHEM_TYPE": "'OBSERVER'",
            "CENTER": f"'{stn}@399'", "TLIST": tlist, "TLIST_TYPE": "'CAL'",
            "TIME_TYPE": "'UT'", "QUANTITIES": "'1'", "REF_SYSTEM": "'ICRF'",
            "ANG_FORMAT": "'DEG'", "CSV_FORMAT": "'YES'", "EXTRA_PREC": "'YES'",
            "CAL_FORMAT": "'BOTH'",
        }
        r = requests.get(HORIZONS_URL, params=params, timeout=90)
        r.raise_for_status()
        result = r.json().get("result", "")
        if "** Error" in result:
            raise RuntimeError(result[:700])
        z = parse_observer(result)
        if len(z) != len(rr):
            raise RuntimeError(f"{stn}: returned {len(z)} for {len(rr)}")
        out.extend(z)
    return out


def oc_rows(name, command):
    obs = fetch_mpc(name)
    groups = defaultdict(list)
    for x in obs:
        groups[x["stn"]].append(x)
    rows, failures = [], []
    for stn in sorted(groups):
        rr = sorted(groups[stn], key=lambda x: x["time"])
        try:
            pp = observer_query(command, stn, rr)
            for o, (jd, rap, decp) in zip(rr, pp, strict=True):
                dra = ((o["ra_obs_deg"] - rap + 180.0) % 360.0) - 180.0
                dra_as = dra * 3600.0 * math.cos(math.radians(decp))
                ddec_as = (o["dec_obs_deg"] - decp) * 3600.0
                rows.append({
                    "object": name,
                    "obs_id": o["obs_id"],
                    "time": o["time"],
                    "jd_ut": jd,
                    "stn": stn,
                    "catalog": o["catalog"],
                    "mode": o["mode"],
                    "ra_obs_deg": o["ra_obs_deg"],
                    "dec_obs_deg": o["dec_obs_deg"],
                    "ra_pred_deg": rap,
                    "dec_pred_deg": decp,
                    "dra_cosdec_arcsec": dra_as,
                    "ddec_arcsec": ddec_as,
                    "rho_arcsec": math.hypot(dra_as, ddec_as),
                })
        except Exception as exc:
            failures.append({"stn": stn, "n": len(rr), "error": type(exc).__name__ + ":" + str(exc)})
    rows.sort(key=lambda x: (x["jd_ut"], x["stn"], str(x["obs_id"])))
    return rows, failures


def nightly(rows):
    groups = defaultdict(list)
    for r in rows:
        groups[(r["stn"], str(r["catalog"] or "UNKNOWN"), r["time"][:10])].append(r)
    out = []
    for (stn, cat, date), rr in groups.items():
        out.append({
            "stn": stn, "catalog": cat, "date": date, "n_raw": len(rr),
            "jd_ut": sum(x["jd_ut"] for x in rr) / len(rr),
            "dra_arcsec": sum(x["dra_cosdec_arcsec"] for x in rr) / len(rr),
            "ddec_arcsec": sum(x["ddec_arcsec"] for x in rr) / len(rr),
        })
    out.sort(key=lambda x: x["jd_ut"])
    return out


def parse_vectors(result):
    if "$$SOE" not in result or "$$EOE" not in result:
        raise RuntimeError("vector response missing SOE/EOE")
    body = result.split("$$SOE", 1)[1].split("$$EOE", 1)[0].strip()
    out = []
    for rec in csv.reader(io.StringIO(body)):
        v = [x.strip() for x in rec]
        if len(v) < 8:
            continue
        q = [finite(v[i]) for i in (0, 2, 3, 4, 5, 6, 7)]
        if any(x is None for x in q):
            raise RuntimeError(f"bad vector row: {v!r}")
        out.append({
            "jd": q[0],
            "r": np.array(q[1:4], dtype=float),
            "v_yr": np.array(q[4:7], dtype=float) * DAY_PER_YEAR,
        })
    return out


def vector_query(command, center, jds):
    out = []
    for i in range(0, len(jds), CHUNK):
        chunk = jds[i:i+CHUNK]
        tlist = " ".join(f"'{x:.12f}'" for x in chunk)
        params = {
            "format": "json", "COMMAND": command, "OBJ_DATA": "'NO'",
            "MAKE_EPHEM": "'YES'", "EPHEM_TYPE": "'VECTORS'",
            "CENTER": center, "TLIST": tlist, "TLIST_TYPE": "'JD'",
            "TIME_TYPE": "'UT'", "REF_SYSTEM": "'ICRF'", "REF_PLANE": "'FRAME'",
            "OUT_UNITS": "'AU-D'", "VEC_TABLE": "'2'", "CSV_FORMAT": "'YES'",
        }
        r = requests.get(HORIZONS_URL, params=params, timeout=90)
        r.raise_for_status()
        result = r.json().get("result", "")
        if "** Error" in result:
            raise RuntimeError(result[:700])
        out.extend(parse_vectors(result))
    if len(out) != len(jds):
        raise RuntimeError(f"vector query returned {len(out)} for {len(jds)}")
    return out


def vec5_to_T(t):
    a,b,c,d,e = map(float, t)
    return np.array([[a,b,c],[b,d,e],[c,e,-a-d]], dtype=float)


def initial_states(command, tref):
    states = {}
    for name, cmd, mass in BODIES:
        z = vector_query(cmd, "'500@10'", [tref])[0]
        states[name] = {"r": z["r"], "v": z["v_yr"], "m": mass}
    z = vector_query(command, "'500@10'", [tref])[0]
    states["TARGET"] = {"r": z["r"], "v": z["v_yr"], "m": 0.0}
    return states


def build_sim(states, target_delta=None, tide_vec=None):
    sim = rebound.Simulation()
    sim.units = ("AU", "yr", "Msun")
    sim.integrator = "ias15"
    sim.ri_ias15.epsilon = 1e-10
    sim.add(m=1.0, hash="Sun")
    for name, _, _ in BODIES:
        s = states[name]
        sim.add(
            m=s["m"], x=s["r"][0], y=s["r"][1], z=s["r"][2],
            vx=s["v"][0], vy=s["v"][1], vz=s["v"][2], hash=name,
        )
    s = states["TARGET"]
    r = s["r"].copy()
    v = s["v"].copy()
    if target_delta is not None:
        r += np.asarray(target_delta[:3], dtype=float)
        v += np.asarray(target_delta[3:], dtype=float)
    sim.add(m=0.0, x=r[0], y=r[1], z=r[2], vx=v[0], vy=v[1], vz=v[2], hash="TARGET")
    sim.move_to_com()
    if tide_vec is not None:
        T = vec5_to_T(tide_vec)
        def tide_force(reb_sim):
            rs = reb_sim.contents
            ps = rs.particles
            sun = ps[0]
            for ii in range(1, rs.N):
                p = ps[ii]
                dx = p.x - sun.x
                dy = p.y - sun.y
                dz = p.z - sun.z
                p.ax += T[0,0]*dx + T[0,1]*dy + T[0,2]*dz
                p.ay += T[1,0]*dx + T[1,1]*dy + T[1,2]*dz
                p.az += T[2,0]*dx + T[2,1]*dy + T[2,2]*dz
        sim.additional_forces = tide_force
        sim._pbh_tide_force = tide_force
    return sim


def rel_now(sim):
    t = sim.particles["TARGET"]
    e = sim.particles["Earth"]
    return np.array([t.x-e.x, t.y-e.y, t.z-e.z], dtype=float)


def sample_rel(states, offsets, target_delta=None, tide_vec=None):
    out = [None] * len(offsets)
    neg = sorted([(x,i) for i,x in enumerate(offsets) if x < 0.0], reverse=True)
    pos = sorted([(x,i) for i,x in enumerate(offsets) if x >= 0.0])
    for seq in (neg, pos):
        if not seq:
            continue
        sim = build_sim(states, target_delta=target_delta, tide_vec=tide_vec)
        for t, i in seq:
            sim.integrate(float(t), exact_finish_time=1)
            out[i] = rel_now(sim)
    return out


def tangent_operator(r):
    rho = float(np.linalg.norm(r))
    u = r / rho
    ra = math.atan2(u[1], u[0])
    dec = math.asin(u[2])
    ea = np.array([-math.sin(ra), math.cos(ra), 0.0])
    ed = np.array([-math.cos(ra)*math.sin(dec), -math.sin(ra)*math.sin(dec), math.cos(dec)])
    return np.vstack([ea, ed]) / rho


def tangent_difference(r_ref, r_other):
    rho = float(np.linalg.norm(r_ref))
    u = r_ref / rho
    v = r_other / float(np.linalg.norm(r_other))
    ra = math.atan2(u[1], u[0])
    dec = math.asin(u[2])
    ea = np.array([-math.sin(ra), math.cos(ra), 0.0])
    ed = np.array([-math.cos(ra)*math.sin(dec), -math.sin(ra)*math.sin(dec), math.cos(dec)])
    return np.array([ea @ (v-u), ed @ (v-u)]) * ARCSEC_PER_RAD


def nuisance_matrix(nights, A):
    n = len(nights)
    blocks = [A]
    cra = np.zeros(2*n); cdec = np.zeros(2*n)
    cra[0::2] = 1.0; cdec[1::2] = 1.0
    blocks += [cra[:,None], cdec[:,None]]
    for key in ("stn", "catalog"):
        levels = sorted({x[key] for x in nights})
        for lev in levels[1:]:
            ra = np.zeros(2*n); de = np.zeros(2*n)
            for i,x in enumerate(nights):
                if x[key] == lev:
                    ra[2*i] = 1.0
                    de[2*i+1] = 1.0
            blocks += [ra[:,None], de[:,None]]
    return np.column_stack(blocks)


def profile_design(nights, A, B, y):
    X = nuisance_matrix(nights, A)
    norms = np.linalg.norm(X, axis=0)
    keep = norms > 1e-15
    Xn = X[:, keep] / norms[keep]
    beta = np.linalg.lstsq(Xn, y, rcond=RCOND)[0]
    gamma = np.linalg.lstsq(Xn, B, rcond=RCOND)[0]
    yp = y - Xn @ beta
    Bp = B - Xn @ gamma
    return Xn, yp, Bp


def object_design(name, command, raw, nights):
    jds = [x["jd_ut"] for x in nights]
    tref = sum(jds)/len(jds)
    offsets = [(jd-tref)/DAY_PER_YEAR for jd in jds]
    states = initial_states(command, tref)
    geo = vector_query(command, "'500@399'", jds)
    M = [tangent_operator(g["r"]) for g in geo]

    print("PBH_NBODY_PROGRESS="+json.dumps({"object":name,"stage":"baseline","n_nights":len(nights)}), flush=True)
    base = sample_rel(states, offsets)
    mism = np.vstack([tangent_difference(g["r"], b) for g,b in zip(geo,base,strict=True)])

    A = np.zeros((2*len(nights), 6))
    for k in range(6):
        delta = np.zeros(6)
        eps = DPOS_AU if k < 3 else DVEL_AU_YR
        delta[k] = eps
        print("PBH_NBODY_PROGRESS="+json.dumps({"object":name,"stage":"A","column":k}), flush=True)
        plus = sample_rel(states, offsets, target_delta=delta)
        minus = sample_rel(states, offsets, target_delta=-delta)
        for i,(rp,rm) in enumerate(zip(plus,minus,strict=True)):
            A[2*i:2*i+2,k] = M[i] @ ((rp-rm)/(2.0*eps))

    B = np.zeros((2*len(nights), 5))
    for k in range(5):
        tv = np.zeros(5); tv[k] = DTIDE_YR2
        print("PBH_NBODY_PROGRESS="+json.dumps({"object":name,"stage":"B","column":k}), flush=True)
        plus = sample_rel(states, offsets, tide_vec=tv)
        minus = sample_rel(states, offsets, tide_vec=-tv)
        for i,(rp,rm) in enumerate(zip(plus,minus,strict=True)):
            B[2*i:2*i+2,k] = M[i] @ ((rp-rm)/(2.0*DTIDE_YR2))

    # One finite-difference step-size consistency gate.
    tv = np.zeros(5); tv[0] = 0.5*DTIDE_YR2
    plus = sample_rel(states, offsets, tide_vec=tv)
    minus = sample_rel(states, offsets, tide_vec=-tv)
    b_half = np.zeros(2*len(nights))
    for i,(rp,rm) in enumerate(zip(plus,minus,strict=True)):
        b_half[2*i:2*i+2] = M[i] @ ((rp-rm)/DTIDE_YR2)
    fd_rel = float(np.linalg.norm(b_half-B[:,0]) / max(np.linalg.norm(B[:,0]),1e-300))

    y = np.zeros(2*len(nights))
    for i,n in enumerate(nights):
        y[2*i] = n["dra_arcsec"]/ARCSEC_PER_RAD
        y[2*i+1] = n["ddec_arcsec"]/ARCSEC_PER_RAD

    Xn, yp, Bp = profile_design(nights, A, B, y)
    s = np.linalg.svd(Bp, compute_uv=False)
    frac = []
    for k in range(5):
        den = float(B[:,k] @ B[:,k])
        frac.append(float(Bp[:,k] @ Bp[:,k])/den if den > 0 else 0.0)
    return {
        "A": A, "B": B, "Bp": Bp, "yp": yp,
        "diag": {
            "n_raw": len(raw), "n_nights": len(nights), "t_ref_jd_ut": tref,
            "rank_nuisance": int(np.linalg.matrix_rank(Xn)),
            "Bp_rank": int(np.linalg.matrix_rank(Bp)),
            "Bp_singular_values": s.tolist(),
            "Bp_condition": float(s[0]/s[-1]) if s[-1] > 0 else float("inf"),
            "tide_column_information_fraction_after_nuisance": frac,
            "tide_fd_q0_step_halving_relative_difference": fd_rel,
            "baseline_vs_horizons_tangent_rms_arcsec": float(np.sqrt(np.mean(mism*mism))),
            "baseline_vs_horizons_tangent_max_abs_arcsec": float(np.max(np.abs(mism))),
            "y_rms_arcsec": float(np.sqrt(np.mean(y*y))*ARCSEC_PER_RAD),
            "yp_rms_arcsec": float(np.sqrt(np.mean(yp*yp))*ARCSEC_PER_RAD),
        },
    }


def point_mass_diag(T):
    vals, vecs = np.linalg.eigh(T)
    axis = vecs[:, -1]
    Q = 3*np.outer(axis,axis)-np.eye(3)
    mu = float(np.sum(T*Q)/6.0)
    res = float(np.linalg.norm(T-mu*Q, "fro")/max(np.linalg.norm(T,"fro"),1e-300))
    split = float(abs(vals[1]-vals[0])/max(abs(mu),1e-300))
    ra = math.degrees(math.atan2(axis[1],axis[0])) % 360.0
    dec = math.degrees(math.asin(axis[2]/np.linalg.norm(axis)))
    return vals, axis, ra, dec, mu, res, split


def main():
    raw_all, raw_by, failures = [], {}, {}
    for name, command in OBJECTS.items():
        print("PBH_NBODY_PROGRESS="+json.dumps({"object":name,"stage":"oc"}), flush=True)
        rows, fail = oc_rows(name, command)
        raw_by[name] = rows; raw_all.extend(rows); failures[name] = fail
    sha = hashlib.sha256(json.dumps(raw_all,sort_keys=True,separators=(",",":")).encode()).hexdigest()
    if sha != EXPECTED_RAW_SHA:
        raise RuntimeError(f"frozen O-C hash mismatch {sha}")
    if any(failures.values()):
        raise RuntimeError("station failures "+json.dumps(failures,sort_keys=True))

    designs = {}
    Bs, ys = [], []
    nightly_payload = {}
    for name, command in OBJECTS.items():
        ns = nightly(raw_by[name]); nightly_payload[name] = ns
        d = object_design(name, command, raw_by[name], ns)
        designs[name] = d
        Bs.append(d["Bp"]); ys.append(d["yp"])

    B = np.vstack(Bs); y = np.concatenate(ys)
    s = np.linalg.svd(B, compute_uv=False)
    t = np.linalg.lstsq(B, y, rcond=RCOND)[0]
    T = vec5_to_T(t)
    pred = B @ t
    rss0 = float(y@y); rss = float((y-pred)@(y-pred))
    vals, axis, ra, dec, mu, shape, split = point_mass_diag(T)
    night_sha = hashlib.sha256(json.dumps(nightly_payload,sort_keys=True,separators=(",",":")).encode()).hexdigest()

    out = {
        "schema": "PBH_MAJOR_PLANET_REBOUND_VARIATIONAL_TIDAL_V0_2",
        "status": "ANALYZED_DIAGNOSTIC_ONLY",
        "software": {"rebound": rebound.__version__, "integrator": "IAS15", "ias15_epsilon": 1e-10},
        "raw_rows_sha256": sha, "nightly_rows_sha256": night_sha,
        "finite_difference": {"dpos_au":DPOS_AU,"dvel_au_yr":DVEL_AU_YR,"dtide_yr2":DTIDE_YR2},
        "objects": {k:v["diag"] for k,v in designs.items()},
        "common": {
            "n_rows": int(B.shape[0]), "rank": int(np.linalg.matrix_rank(B)),
            "singular_values": s.tolist(),
            "condition": float(s[0]/s[-1]) if s[-1] > 0 else float("inf"),
            "tensor_vec5_yr2": t.tolist(), "tensor_3x3_yr2": T.tolist(),
            "eigenvalues_yr2": vals.tolist(),
            "axis_largest_eigenvalue": {"ra_deg":ra,"dec_deg":dec,"vector":axis.tolist()},
            "mu_frobenius_yr2": mu,
            "point_mass_shape_residual": shape,
            "transverse_split_over_mu": split,
            "nuisance_profiled_rms_arcsec": float(math.sqrt(rss0/len(y))*ARCSEC_PER_RAD),
            "post_tide_rms_arcsec": float(math.sqrt(rss/len(y))*ARCSEC_PER_RAD),
            "explained_fraction": 1.0-rss/rss0 if rss0 > 0 else 0.0,
        },
        "model": {
            "massive_bodies": [x[0] for x in BODIES],
            "external_tide_applied_relative_to_sun_to_all_non_sun_particles": True,
            "target_massless": True,
            "observer_dynamics": "Earth center in same N-body+tide simulation",
        },
        "interpretation_boundary": {
            "planet_masses_approximate": True,
            "no_minor_planet_perturbers": True,
            "no_relativistic_force": True,
            "no_light_time_derivative": True,
            "topocentric_station_displacement_not_in_derivative": True,
            "nightly_equal_weight": True,
            "current_JPL_not_independent_of_MPC": True,
            "PBH_not_inferred": True,
        },
        "canon_allowed": False,
    }
    print("PBH_NBODY_SUMMARY="+json.dumps(out,sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
