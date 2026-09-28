#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import math
from collections import Counter, defaultdict
from datetime import datetime, timedelta, timezone

import numpy as np
import requests

CATALOG_URL = (
    "https://vizier.cds.unistra.fr/viz-bin/asu-tsv"
    "?-source=J%2FAJ%2F152%2F80%2Fobs"
    "&-out.all&-out.max=unlimited"
)
PLUTO_PACKED = "D4340"
PLUTO_PERMID = "134340"


def parse_catalog(text: str) -> list[dict[str,str]]:
    lines=[line for line in text.splitlines() if line and not line.startswith("#")]
    header_idx=None
    for i,line in enumerate(lines):
        if "Planet" in line and "Date" in line and "RAJ2000" in line and "DEJ2000" in line:
            header_idx=i; break
    if header_idx is None:
        raise RuntimeError("catalog header not found")
    header=lines[header_idx].split("\t")
    rows=[]
    for line in lines[header_idx+1:]:
        parts=line.split("\t")
        if len(parts)!=len(header): continue
        row=dict(zip(header,parts))
        if not row.get("Date","").strip()[:4].isdigit(): continue
        rows.append(row)
    return rows


def mpc_date_to_iso(text: str) -> str:
    y,m,d=text.strip().split(":")
    year=int(y); month=int(m); dayf=float(d)
    day=int(math.floor(dayf))
    frac=dayf-day
    base=datetime(year,month,day,tzinfo=timezone.utc)+timedelta(days=frac)
    # ADES permits fractional seconds. Microseconds preserve the source precision.
    return base.strftime("%Y-%m-%dT%H:%M:%S.") + f"{base.microsecond:06d}Z"


def ra_hms_to_deg(text: str) -> float:
    h,m,s=text.strip().split()
    return 15.0*(float(h)+float(m)/60.0+float(s)/3600.0)


def dec_dms_to_deg(text: str) -> float:
    parts=text.strip().split()
    sign=-1.0 if parts[0].startswith("-") else 1.0
    d=abs(float(parts[0])); m=float(parts[1]); s=float(parts[2])
    return sign*(d+m/60.0+s/3600.0)


def mode_from_legacy(code: str) -> str:
    code=code.strip()
    return {
        "C":"CCD", "c":"CCD",
        "E":"OCC",
        "e":"ENC",
    }.get(code,"PHO")


def make_psv(rows: list[dict[str,str]]) -> str:
    header="permID|mode|stn|obsTime|ra|dec|rmsRA|rmsDec|mag|band"
    out=["# version=2022",header]
    for row in rows:
        mag=row.get("mag","").strip()
        band=row.get("Flt","").strip()
        rmsra=row.get("e_RAJ2000","").strip()
        rmsdec=row.get("e_DEJ2000","").strip()
        vals=[
            PLUTO_PERMID,
            mode_from_legacy(row.get("R","")),
            row.get("Obs","").strip(),
            mpc_date_to_iso(row["Date"]),
            f"{ra_hms_to_deg(row['RAJ2000']):.10f}",
            f"{dec_dms_to_deg(row['DEJ2000']):.10f}",
            rmsra,
            rmsdec,
            mag,
            band,
        ]
        out.append("|".join(vals))
    return "\n".join(out)+"\n"


def stable_rows_sha(rows: list[dict[str,str]]) -> str:
    fields=["Planet","Prov","Date","RAJ2000","DEJ2000","Obs","e_RAJ2000","e_DEJ2000","Out"]
    payload="\n".join("\t".join(row.get(k,"").strip() for k in fields) for row in rows).encode()
    return hashlib.sha256(payload).hexdigest()


def arr_col(table, name: str, dtype=float):
    col=getattr(table,name)
    try:
        return np.asarray(col.to_numpy(zero_copy_only=False),dtype=dtype)
    except TypeError:
        return np.asarray(col.to_numpy(),dtype=dtype)


def list_col(table,name: str):
    col=getattr(table,name)
    return col.to_pylist()


def bin_residuals(year, ra, dec, selected):
    edges=[1930,1935,1940,1945,1950,1960,1980,2000,2010,2016]
    out=[]
    for lo,hi in zip(edges[:-1],edges[1:]):
        m=(year>=lo)&(year<hi)&selected
        if not np.any(m):
            continue
        out.append({
            "year_range":[lo,hi],
            "n":int(np.sum(m)),
            "ra_mean_arcsec":float(np.mean(ra[m])),
            "ra_median_arcsec":float(np.median(ra[m])),
            "ra_rms_arcsec":float(np.sqrt(np.mean(ra[m]**2))),
            "dec_mean_arcsec":float(np.mean(dec[m])),
            "dec_median_arcsec":float(np.median(dec[m])),
            "dec_rms_arcsec":float(np.sqrt(np.mean(dec[m]**2))),
        })
    return out


def slope_report(year, values, selected, lo, hi):
    m=(year>=lo)&(year<hi)&selected&np.isfinite(values)
    if np.sum(m)<3:
        return {"year_range":[lo,hi],"n":int(np.sum(m)),"slope_arcsec_per_year":None}
    x=year[m]
    y=values[m]
    x0=float(np.mean(x))
    slope,intercept=np.polyfit(x-x0,y,1)
    pred=slope*(x-x0)+intercept
    rms=float(np.sqrt(np.mean((y-pred)**2)))
    return {
        "year_range":[lo,hi],
        "n":int(np.sum(m)),
        "slope_arcsec_per_year":float(slope),
        "intercept_at_mean_epoch_arcsec":float(intercept),
        "mean_epoch":x0,
        "detrended_rms_arcsec":rms,
    }


def main() -> int:
    print("PLUTO_OD_STAGE=FETCH_CDS",flush=True)
    response=requests.get(
        CATALOG_URL,timeout=180,
        headers={"User-Agent":"GREMLIN-FERRET-PBH-ORBITAL-OD/0.1 read-only"}
    )
    response.raise_for_status()
    rows=parse_catalog(response.text)
    pluto_all=[r for r in rows if r.get("Planet","").strip()==PLUTO_PACKED]
    pluto=[r for r in pluto_all if r.get("Out","").strip()!="1"]
    if len(pluto_all)!=6677:
        raise RuntimeError(f"expected 6677 Pluto rows, got {len(pluto_all)}")
    print(f"PLUTO_OD_STAGE=CATALOG_OK all={len(pluto_all)} selected={len(pluto)}",flush=True)

    psv=make_psv(pluto)
    print(f"PLUTO_OD_STAGE=PSV_READY bytes={len(psv.encode())}",flush=True)

    import empyrean
    from empyrean import ODConfig, SolveForParams

    print("PLUTO_OD_STAGE=EMPYREAN_INIT",flush=True)
    empyrean.initialize()
    print("PLUTO_OD_STAGE=EMPYREAN_READY version="+empyrean.version_string(),flush=True)

    obs,radar=empyrean.read_ades(psv)
    if len(obs)!=len(pluto):
        raise RuntimeError(f"ADES row count mismatch {len(obs)} != {len(pluto)}")

    cfg=ODConfig(
        solve_for=SolveForParams.STATE_ONLY,
        allow_arc_truncation=False,
        fit_station_biases=False,
        num_threads=0,
    )
    print("PLUTO_OD_STAGE=FIT_START",flush=True)
    fits=empyrean.determine(obs,config=cfg)
    fit=fits.single()
    print(
        f"PLUTO_OD_STAGE=FIT_DONE converged={fit.converged} "
        f"termination={fit.termination} chi2={fit.summary.reduced_chi2}",
        flush=True,
    )

    res=fit.observations
    epoch=arr_col(res,"epoch_mjd_tdb")
    ra=arr_col(res,"ra_residual")
    dec=arr_col(res,"dec_residual")
    selected=np.asarray(list_col(res,"selected"),dtype=bool)
    codes=list_col(res,"obs_code")
    year=2000.0+(epoch-51544.5)/365.25

    station=defaultdict(lambda:{"n":0,"ra":[],"dec":[],"years":[]})
    for i,code in enumerate(codes):
        if not selected[i]: continue
        row=station[str(code)]
        row["n"]+=1; row["ra"].append(float(ra[i])); row["dec"].append(float(dec[i])); row["years"].append(float(year[i]))
    station_summary=[]
    for code,row in station.items():
        station_summary.append({
            "obs_code":code,
            "n":row["n"],
            "year_min":min(row["years"]),
            "year_max":max(row["years"]),
            "ra_mean_arcsec":float(np.mean(row["ra"])),
            "dec_mean_arcsec":float(np.mean(row["dec"])),
            "ra_rms_arcsec":float(np.sqrt(np.mean(np.asarray(row["ra"])**2))),
            "dec_rms_arcsec":float(np.sqrt(np.mean(np.asarray(row["dec"])**2))),
        })
    station_summary.sort(key=lambda x:x["n"],reverse=True)

    orbit=fit.orbit
    orbit_dict={"representation":type(orbit).__name__}
    # Keep the delivered state/covariance provenance without assuming a column flavor.
    try:
        orbit_dict["table"]=orbit.table.to_pylist()
    except Exception:
        try: orbit_dict["repr"]=repr(orbit)
        except Exception: pass

    receipt={
        "schema":"PBH_PLUTO_BASELINE_OD_V0_1",
        "status":"FITTED_BASELINE_NO_EXTERNAL_PERTURBER",
        "source":{
            "catalog":"CDS/VizieR J/AJ/152/80/obs",
            "source_url":CATALOG_URL,
            "raw_response_sha256":hashlib.sha256(response.content).hexdigest(),
            "normalized_pluto_rows_sha256":stable_rows_sha(pluto_all),
            "pluto_packed_id":PLUTO_PACKED,
            "pluto_perm_id":PLUTO_PERMID,
            "n_pluto_catalog":len(pluto_all),
            "n_preflagged_outliers":len(pluto_all)-len(pluto),
            "n_input_fit":len(pluto),
        },
        "engine":{
            "name":"empyrean",
            "version":empyrean.version_string(),
            "force_model":str(fit.force_model_used),
            "solve_for":str(fit.solve_for_used),
            "converged":bool(fit.converged),
            "termination":fit.termination,
            "iterations":int(fit.iterations),
            "final_solve_iterations":fit.final_solve_iterations,
            "accepted_steps":int(fit.accepted_steps),
            "warnings":list(fit.warnings),
        },
        "fit_summary":{
            "reduced_chi2":float(fit.summary.reduced_chi2),
            "rms_ra_arcsec":float(fit.summary.rms_ra_arcsec),
            "rms_dec_arcsec":float(fit.summary.rms_dec_arcsec),
            "num_obs":int(fit.summary.num_obs),
            "num_selected":int(fit.summary.num_selected),
            "fit_acceptable":bool(fit.acceptability.fit_acceptable),
            "extrapolation_acceptable":bool(fit.acceptability.extrapolation_acceptable),
        },
        "residual_bins":bin_residuals(year,ra,dec,selected),
        "early_declination_slope":slope_report(year,dec,selected,1930,1952),
        "early_ra_slope":slope_report(year,ra,selected,1930,1952),
        "full_declination_slope":slope_report(year,dec,selected,float(np.min(year)),float(np.max(year))+1e-6),
        "stations":station_summary[:30],
        "orbit":orbit_dict,
        "inference_boundary":{
            "external_tidal_tensor_fitted":False,
            "PBH_fitted":False,
            "meaning":"This run establishes the no-external-perturber residual field only."
        },
        "authority":"OBSERVATIONAL_BASELINE_ONLY",
        "canon_allowed":False,
    }
    compact={
        "schema":receipt["schema"],
        "status":receipt["status"],
        "source":receipt["source"],
        "engine":receipt["engine"],
        "fit_summary":receipt["fit_summary"],
        "early_declination_slope":receipt["early_declination_slope"],
        "early_ra_slope":receipt["early_ra_slope"],
        "residual_bins":receipt["residual_bins"],
        "stations":receipt["stations"][:10],
    }
    print("PBH_PLUTO_OD_SUMMARY="+json.dumps(compact,sort_keys=True,separators=(",",":")),flush=True)
    print("PBH_PLUTO_OD_RECEIPT="+json.dumps(receipt,sort_keys=True,separators=(",",":")),flush=True)
    return 0


if __name__=="__main__":
    raise SystemExit(main())
