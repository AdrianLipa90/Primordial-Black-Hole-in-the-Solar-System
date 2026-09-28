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

CHUNK = 20
ARCSEC_TO_RAD = math.pi / (180.0 * 3600.0)
SECONDS_PER_YEAR = 365.25 * 86400.0
GAIA_CATALOGS = {"Gaia1", "Gaia2", "Gaia3", "Gaia3E"}


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


def horizons_batch(command: str, station: str, rows: list[dict]):
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
    last_error = None
    for attempt in range(4):
        try:
            response = requests.get(HORIZONS_URL, params=params, timeout=90)
            if response.status_code in (429, 502, 503, 504):
                raise requests.HTTPError(
                    f"{response.status_code} transient Horizons response",
                    response=response,
                )
            response.raise_for_status()
            payload = response.json()
            break
        except Exception as exc:
            last_error = exc
            if attempt == 3:
                raise
            time.sleep(0.6 * (attempt + 1))
    else:
        raise last_error

    result = payload.get("result", "")
    if "$$SOE" not in result or "$$EOE" not in result:
        raise RuntimeError("Horizons response lacks ephemeris markers: " + result[:500])

    block = result.split("$$SOE", 1)[1].split("$$EOE", 1)[0]
    pred = []
    for line in block.splitlines():
        line = line.strip()
        if not line:
            continue
        nums = []
        for part in (p.strip() for p in line.split(",")):
            try:
                nums.append(float(part))
            except Exception:
                pass
        if len(nums) >= 3:
            pred.append({"jd": nums[0], "ra_deg": nums[1], "dec_deg": nums[2]})

    if len(pred) != len(rows):
        raise RuntimeError(
            f"Horizons row mismatch station={station}: expected {len(rows)}, got {len(pred)}"
        )

    return pred


def fetch_residuals(spec):
    response = requests.get(
        MPC_URL,
        json={"desigs": [spec["designation"]], "output_format": ["ADES_DF"], "ades_version": "2022"},
        timeout=90,
    )
    response.raise_for_status()
    raw_rows = response.json()[0]["ADES_DF"]

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
            "astcat": row.get("astcat") or "",
            "mode": row.get("mode"),
            "ra_obs_deg": ra,
            "dec_obs_deg": dec,
        })

    by_station = collections.defaultdict(list)
    for row in usable:
        by_station[row["stn"]].append(row)

    residuals = []
    failures = []
    for stn in sorted(by_station):
        rows = sorted(by_station[stn], key=lambda r: r["obstime"])
        for part in chunks(rows, CHUNK):
            try:
                pred = horizons_batch(spec["command"], stn, part)
                for obs, ep in zip(part, pred):
                    dra_deg = wrap_deg(obs["ra_obs_deg"] - ep["ra_deg"])
                    cosd = math.cos(math.radians(ep["dec_deg"]))
                    residuals.append({
                        **obs,
                        "jd": ep["jd"],
                        "ra_pred_deg": ep["ra_deg"],
                        "dec_pred_deg": ep["dec_deg"],
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
            time.sleep(0.05)

    residuals.sort(key=lambda r: r["jd"])
    return raw_rows, usable, residuals, failures, dict(skipped)


def solve_2x2(a, b, c, u, v):
    det = a*c - b*b
    if abs(det) <= 1e-30:
        raise ValueError("singular 2x2 normal matrix")
    return ((u*c - v*b)/det, (v*a - u*b)/det), det


def fixed_effect_quadratic(rows, key):
    if len(rows) < 8:
        raise ValueError("too few rows for curvature fit")

    tref = sum(r["jd"] for r in rows) / len(rows)
    raw = []
    groups = collections.defaultdict(list)
    for i, r in enumerate(rows):
        tau = (r["jd"] - tref) / 365.25
        x1 = tau
        x2 = 0.5 * tau * tau
        y = float(r[key])
        g = (r["stn"], r.get("astcat") or "")
        raw.append((x1, x2, y, g))
        groups[g].append(i)

    gp_means = {}
    for g, ids in groups.items():
        gp_means[g] = (
            sum(raw[i][0] for i in ids)/len(ids),
            sum(raw[i][1] for i in ids)/len(ids),
            sum(raw[i][2] for i in ids)/len(ids),
        )

    z = []
    for x1, x2, y, g in raw:
        m1, m2, my = gp_means[g]
        z.append((x1-m1, x2-m2, y-my))

    s11 = sum(x1*x1 for x1,x2,y in z)
    s12 = sum(x1*x2 for x1,x2,y in z)
    s22 = sum(x2*x2 for x1,x2,y in z)
    q1 = sum(x1*y for x1,x2,y in z)
    q2 = sum(x2*y for x1,x2,y in z)
    (b1, b2), det = solve_2x2(s11,s12,s22,q1,q2)

    resid = [y-b1*x1-b2*x2 for x1,x2,y in z]
    dof = max(1, len(z)-len(groups)-2)
    sigma2 = sum(e*e for e in resid)/dof
    inv11 = s22/det
    inv22 = s11/det
    inv12 = -s12/det
    se1 = math.sqrt(max(0.0, sigma2*inv11))
    se2 = math.sqrt(max(0.0, sigma2*inv22))

    return {
        "tref_jd": tref,
        "n": len(rows),
        "n_station_catalog_groups": len(groups),
        "linear_arcsec_per_year": b1,
        "linear_se": se1,
        "curvature_arcsec_per_year2": b2,
        "curvature_se": se2,
        "residual_rms_arcsec": math.sqrt(sum(e*e for e in resid)/len(resid)),
    }


def leave_one_group_out_curvature(rows, key):
    groups = sorted(set((r["stn"], r.get("astcat") or "") for r in rows))
    vals = []
    failures = 0
    for g in groups:
        subset = [r for r in rows if (r["stn"], r.get("astcat") or "") != g]
        try:
            vals.append(fixed_effect_quadratic(subset, key)["curvature_arcsec_per_year2"])
        except Exception:
            failures += 1
    if not vals:
        return {"n": 0, "failures": failures}
    med = statistics.median(vals)
    mad = statistics.median(abs(x-med) for x in vals)
    return {
        "n": len(vals),
        "failures": failures,
        "median": med,
        "mad": mad,
        "min": min(vals),
        "max": max(vals),
        "sign_fraction_positive": sum(x>0 for x in vals)/len(vals),
    }


def unit_sky_from_rows(rows):
    x=y=z=0.0
    for r in rows:
        a=math.radians(r["ra_pred_deg"])
        d=math.radians(r["dec_pred_deg"])
        cd=math.cos(d)
        x += cd*math.cos(a)
        y += cd*math.sin(a)
        z += math.sin(d)
    norm=math.sqrt(x*x+y*y+z*z)
    n=(x/norm,y/norm,z/norm)
    ra=math.degrees(math.atan2(n[1],n[0]))%360.0
    dec=math.degrees(math.asin(max(-1.0,min(1.0,n[2]))))
    ea=(-math.sin(math.radians(ra)), math.cos(math.radians(ra)), 0.0)
    er=math.radians(ra); dr=math.radians(dec)
    ed=(-math.sin(dr)*math.cos(er), -math.sin(dr)*math.sin(er), math.cos(dr))
    return n, ea, ed, ra, dec


BASIS = [
    ((1,0,0),(0,0,0),(0,0,-1)),
    ((0,1,0),(1,0,0),(0,0,0)),
    ((0,0,1),(0,0,0),(1,0,0)),
    ((0,0,0),(0,1,0),(0,0,-1)),
    ((0,0,0),(0,0,1),(0,1,0)),
]


def mv(M,v):
    return tuple(sum(M[i][j]*v[j] for j in range(3)) for i in range(3))


def dot(a,b):
    return sum(x*y for x,y in zip(a,b))


def design_row(tangent,n):
    return [dot(tangent,mv(B,n)) for B in BASIS]


def gaussian_solve(A,b):
    n=len(b)
    aug=[list(map(float,A[i]))+[float(b[i])] for i in range(n)]
    pivots=[]
    for col in range(n):
        p=max(range(col,n),key=lambda r: abs(aug[r][col]))
        if abs(aug[p][col]) < 1e-24:
            raise ValueError("singular normal system")
        aug[col],aug[p]=aug[p],aug[col]
        piv=aug[col][col]
        pivots.append(abs(piv))
        for j in range(col,n+1):
            aug[col][j]/=piv
        for r in range(n):
            if r==col:
                continue
            f=aug[r][col]
            if f==0:
                continue
            for j in range(col,n+1):
                aug[r][j]-=f*aug[col][j]
    return [aug[i][n] for i in range(n)], {
        "min_pivot": min(pivots),
        "max_pivot": max(pivots),
        "pivot_ratio": max(pivots)/max(min(pivots),1e-300),
    }


def tensor_from_vec(t):
    return (
        (t[0],t[1],t[2]),
        (t[1],t[3],t[4]),
        (t[2],t[4],-t[0]-t[3]),
    )


def frob(A):
    return math.sqrt(sum(A[i][j]*A[i][j] for i in range(3) for j in range(3)))


def sub(A,B):
    return tuple(tuple(A[i][j]-B[i][j] for j in range(3)) for i in range(3))


def scale(A,s):
    return tuple(tuple(s*x for x in row) for row in A)


def jacobi_eigh3(M, iterations=96):
    a=[list(map(float,row)) for row in M]
    v=[[1.0 if i==j else 0.0 for j in range(3)] for i in range(3)]
    for _ in range(iterations):
        opts=[(abs(a[0][1]),0,1),(abs(a[0][2]),0,2),(abs(a[1][2]),1,2)]
        largest,p,q=max(opts)
        scale0=max(1.0,max(abs(a[i][j]) for i in range(3) for j in range(3)))
        if largest <= 1e-15*scale0:
            break
        app,aqq,apq=a[p][p],a[q][q],a[p][q]
        phi=0.5*math.atan2(2*apq,aqq-app)
        c,s=math.cos(phi),math.sin(phi)
        for k in range(3):
            if k in (p,q):
                continue
            akp,akq=a[k][p],a[k][q]
            a[k][p]=a[p][k]=c*akp-s*akq
            a[k][q]=a[q][k]=s*akp+c*akq
        a[p][p]=c*c*app-2*s*c*apq+s*s*aqq
        a[q][q]=s*s*app+2*s*c*apq+c*c*aqq
        a[p][q]=a[q][p]=0.0
        for k in range(3):
            vkp,vkq=v[k][p],v[k][q]
            v[k][p]=c*vkp-s*vkq
            v[k][q]=s*vkp+c*vkq
    vals=[a[i][i] for i in range(3)]
    vecs=[]
    for col in range(3):
        q=tuple(v[row][col] for row in range(3))
        norm=math.sqrt(dot(q,q))
        vecs.append(tuple(x/norm for x in q))
    order=sorted(range(3),key=lambda i: vals[i])
    return [vals[i] for i in order],[vecs[i] for i in order]


def point_mass_diag(T):
    vals,vecs=jacobi_eigh3(T)
    axis=vecs[-1]
    for x in axis:
        if abs(x)>1e-15:
            if x<0:
                axis=tuple(-z for z in axis)
            break
    Q=tuple(tuple(3*axis[i]*axis[j]-(1.0 if i==j else 0.0) for j in range(3)) for i in range(3))
    mu=sum(T[i][j]*Q[i][j] for i in range(3) for j in range(3))/6.0
    point=scale(Q,mu)
    shape=frob(sub(T,point))/max(frob(T),1e-300)
    split=abs(vals[1]-vals[0])/max(abs(mu),1e-300)
    ra=math.degrees(math.atan2(axis[1],axis[0]))%360.0
    dec=math.degrees(math.asin(max(-1.0,min(1.0,axis[2]))))
    return {
        "eigenvalues_per_year2": vals,
        "mu_per_year2": mu,
        "axis_sign_degenerate": list(axis),
        "axis_ra_deg_sign_choice": ra,
        "axis_dec_deg_sign_choice": dec,
        "point_mass_shape_residual": shape,
        "transverse_eigenvalue_split_over_mu": split,
    }


def fit_stf(entries):
    rows=[]
    ys=[]
    ws=[]
    for e in entries:
        for basis_name,c,se in [
            ("alpha",e["curvature_ra_rad_per_year2"],e["curvature_ra_se_rad_per_year2"]),
            ("delta",e["curvature_dec_rad_per_year2"],e["curvature_dec_se_rad_per_year2"]),
        ]:
            tangent=e["e_alpha"] if basis_name=="alpha" else e["e_delta"]
            rows.append(design_row(tangent,e["n"]))
            ys.append(c)
            floor=max(abs(c)*0.05,1e-18)
            ws.append(1.0/max(se*se,floor*floor))

    normal=[[0.0]*5 for _ in range(5)]
    rhs=[0.0]*5
    for x,y,w in zip(rows,ys,ws):
        for i in range(5):
            rhs[i]+=w*x[i]*y
            for j in range(5):
                normal[i][j]+=w*x[i]*x[j]

    t,conditioning=gaussian_solve(normal,rhs)
    T=tensor_from_vec(t)
    pred=[sum(x[i]*t[i] for i in range(5)) for x in rows]
    wrss=sum(w*(y-p)*(y-p) for y,p,w in zip(ys,pred,ws))
    rss=sum((y-p)**2 for y,p in zip(ys,pred))
    diag=point_mass_diag(T)
    sec2=SECONDS_PER_YEAR**2
    return {
        "tensor_vec5_per_year2": t,
        "tensor_3x3_per_year2": [list(r) for r in T],
        "tensor_3x3_per_s2": [[x/sec2 for x in r] for r in T],
        "weighted_rss": wrss,
        "unweighted_rms_rad_per_year2": math.sqrt(rss/len(ys)),
        "normal_pivot_condition_proxy": conditioning,
        "point_mass_projection": {
            **diag,
            "mu_per_s2": diag["mu_per_year2"]/sec2,
        },
        "dof_nominal_without_nuisance_uncertainty": len(ys)-5,
    }


def analyze_subset(residuals, subset_name):
    if subset_name=="ALL":
        rows=list(residuals)
    elif subset_name=="GAIA":
        rows=[r for r in residuals if (r.get("astcat") or "") in GAIA_CATALOGS]
    else:
        raise ValueError(subset_name)
    if len(rows)<8:
        return {"status":"INSUFFICIENT_DATA","n":len(rows)}

    fra=fixed_effect_quadratic(rows,"dra_cosdec_arcsec")
    fde=fixed_effect_quadratic(rows,"ddec_arcsec")
    loo_ra=leave_one_group_out_curvature(rows,"dra_cosdec_arcsec")
    loo_de=leave_one_group_out_curvature(rows,"ddec_arcsec")
    n,ea,ed,ra,dec=unit_sky_from_rows(rows)
    return {
        "status":"FIT",
        "n":len(rows),
        "sky_mean_ra_deg":ra,
        "sky_mean_dec_deg":dec,
        "n_vec":list(n),
        "e_alpha":list(ea),
        "e_delta":list(ed),
        "ra_fit":fra,
        "dec_fit":fde,
        "loo_station_catalog": {"ra":loo_ra,"dec":loo_de},
        "curvature_ra_rad_per_year2":fra["curvature_arcsec_per_year2"]*ARCSEC_TO_RAD,
        "curvature_ra_se_rad_per_year2":fra["curvature_se"]*ARCSEC_TO_RAD,
        "curvature_dec_rad_per_year2":fde["curvature_arcsec_per_year2"]*ARCSEC_TO_RAD,
        "curvature_dec_se_rad_per_year2":fde["curvature_se"]*ARCSEC_TO_RAD,
    }


def main():
    raw = {}
    for spec in OBJECTS:
        raw_rows, usable, residuals, failures, skipped = fetch_residuals(spec)
        canonical=json.dumps(residuals,sort_keys=True,separators=(",",":")).encode()
        packet={
            "designation":spec["designation"],
            "status":"PASS" if residuals and not failures else ("PARTIAL" if residuals else "FAIL"),
            "n_mpc_rows":len(raw_rows),
            "n_usable":len(usable),
            "n_residual_rows":len(residuals),
            "failures":failures,
            "skipped":skipped,
            "residual_sha256":hashlib.sha256(canonical).hexdigest(),
            "ALL":analyze_subset(residuals,"ALL"),
            "GAIA":analyze_subset(residuals,"GAIA"),
        }
        raw[spec["designation"]]=packet
        print("PBH_TIDAL_SCREEN_OBJECT="+json.dumps(packet,sort_keys=True),flush=True)

    modes={}
    for subset in ("ALL","GAIA"):
        entries=[]
        for spec in OBJECTS:
            e=raw[spec["designation"]][subset]
            if e.get("status")=="FIT":
                entries.append(e)
        modes[subset]={
            "n_tracers":len(entries),
            "fit":fit_stf(entries) if len(entries)>=3 else None,
        }

    receipt={
        "schema":"PBH_THREE_TRACER_STF_GEOMETRY_SCREEN_V0_3",
        "status":"SCREENING_ONLY",
        "authority":"OBSERVED_MINUS_JPL_GEOMETRY_SCREEN",
        "objects":raw,
        "modes":modes,
        "method":{
            "residual_model":"per coordinate: stationxcatalog fixed effects + global linear term + 0.5*curvature*t^2",
            "geometry":"theta_ddot ~= tangent^T T n",
            "tensor":"generic symmetric trace-free 3x3, 5 parameters",
            "point_mass_test":"eigenshape projection only",
            "leave_one_group_out":"stationxcatalog curvature robustness diagnostic",
        },
        "interpretation_boundary":{
            "JPL_orbits_fit_same_MPC_data":True,
            "curvature_not_independent_holdout":True,
            "local_6D_orbit_state_not_variationally_profiled":True,
            "three_tracers_give_only_one_nominal_tensor_dof":True,
            "short_arc_and_catalog_systematics_can_mimic_curvature":True,
            "no_sigma_claim":True,
            "no_PBH_identity_claim":True,
            "historical_PBH_coordinate_not_used":True,
            "historical_576_AU_not_used":True,
        },
        "canon_allowed":False,
    }
    print("PBH_TIDAL_SCREEN_RECEIPT="+json.dumps(receipt,sort_keys=True),flush=True)


if __name__=="__main__":
    main()
