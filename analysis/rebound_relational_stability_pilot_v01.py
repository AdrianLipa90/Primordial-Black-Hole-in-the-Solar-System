#!/usr/bin/env python3
from __future__ import annotations

import json, math, urllib.parse, urllib.request
try:
    import rebound
except ModuleNotFoundError:
    import subprocess, sys
    subprocess.check_call([sys.executable, '-m', 'pip', 'install', '--no-cache-dir', 'rebound'])
    import rebound

MU_REL = 2.690873867e-7  # Earth masses / AU^3 from a_crit=170 AU leading-order balance
M_EARTH_PER_MSUN = 332946.0
VARPI_P_DEG = 242.1877460640173
T_YR = 200_000.0
DT_YR = 1.0
PHASES_DEG = (0.0, 120.0, 240.0)

DESIGNATIONS = [
"2016 TP120","2016 QV89","2015 UN105","2003 SS422","2013 UT15",
"2000 CR105","2012 VP113","2013 FT28","2018 VM35","2014 WB556",
"2014 SR349","2010 GB174","2004 VN112","2016 SD106","2015 RX245",
"2015 BP519","2013 RA109","2003 VB12","2021 DK18","2013 SY99",
]

GIANTS = [
    ("Jupiter",9.5458e-4,5.204),
    ("Saturn",2.858e-4,9.583),
    ("Uranus",4.366e-5,19.191),
    ("Neptune",5.151e-5,30.07),
]

A_GRID = (250.0, 290.0, 330.0, 380.0)
E_GRID = (0.10, 0.30, 0.50)

def wrap(x):
    return (x + 180.0) % 360.0 - 180.0

def mass_on_mu(a,e):
    return MU_REL * a**3 * (1.0-e*e)**1.5

def fetch_tnos():
    fields="full_name,epoch,a,e,q,i,om,w,ma,condition_code,data_arc,n_obs_used"
    cdata=json.dumps({"AND":["a|GT|150","q|GT|30"]},separators=(",",":"))
    params={"fields":fields,"sb-class":"TNO","sb-cdata":cdata,"limit":"500","full-prec":"true"}
    url="https://ssd-api.jpl.nasa.gov/sbdb_query.api?"+urllib.parse.urlencode(params)
    req=urllib.request.Request(url,headers={"User-Agent":"PBH-relational-rebound/0.1"})
    with urllib.request.urlopen(req,timeout=60) as resp:
        payload=json.loads(resp.read().decode("utf-8"))
    rows=[]
    names=payload["fields"]
    for raw in payload["data"]:
        d=dict(zip(names,raw))
        full=str(d["full_name"])
        des=next((x for x in DESIGNATIONS if x in full),None)
        if des is None:
            continue
        try:
            row={
                "designation":des,
                "full_name":full.strip(),
                "epoch":float(d["epoch"]),
                "a":float(d["a"]),"e":float(d["e"]),"q":float(d["q"]),
                "i":float(d["i"]),"Omega":float(d["om"]),"omega":float(d["w"]),
                "M":float(d["ma"]),
            }
        except Exception:
            continue
        rows.append(row)
    found={x["designation"] for x in rows}
    missing=[x for x in DESIGNATIONS if x not in found]
    if missing:
        raise RuntimeError("missing JPL rows: "+",".join(missing))
    rows.sort(key=lambda x: DESIGNATIONS.index(x["designation"]))
    return url,rows

def phase_order(sim, start_idx):
    z=0j; wsum=0.0; qs=[]; bound=0
    rows=[]
    sun=sim.particles[0]
    for idx in range(start_idx,len(sim.particles)):
        p=sim.particles[idx]
        try:
            o=p.orbit(primary=sun)
        except Exception:
            continue
        if not math.isfinite(o.a) or o.a<=0 or o.e>=1:
            rows.append({"bound":False})
            continue
        varpi=(math.degrees(o.Omega+o.omega))%360.0
        q=o.a*(1-o.e)
        w=max(1e-6,o.e)
        z += w*complex(math.cos(math.radians(varpi)),math.sin(math.radians(varpi)))
        wsum += w
        qs.append(q); bound+=1
        rows.append({"bound":True,"a":o.a,"e":o.e,"q":q,"varpi":varpi})
    R=abs(z)/wsum if wsum else 0.0
    mu=(math.degrees(math.atan2(z.imag,z.real))%360.0) if wsum else None
    return {"R":R,"mean_varpi_deg":mu,"bound":bound,"min_q_au":min(qs) if qs else None,"rows":rows}

def make_sim(tnos,a_p,e_p,m_earth,M_p_deg):
    sim=rebound.Simulation()
    sim.units=("AU","yr","Msun")
    sim.add(m=1.0)
    # Circular-coplanar giant-planet secular baseline; mean longitudes staggered.
    for k,(name,m,a) in enumerate(GIANTS):
        sim.add(m=m,a=a,e=0.0,inc=0.0,Omega=0.0,omega=0.0,M=math.radians((k*90.0)%360))
    sim.add(
        m=m_earth/M_EARTH_PER_MSUN,
        a=a_p,e=e_p,inc=0.0,Omega=0.0,
        omega=math.radians(VARPI_P_DEG),
        M=math.radians(M_p_deg),
        hash="PBH_LATENT"
    )
    start_idx=len(sim.particles)
    for row in tnos:
        sim.add(
            m=0.0,a=row["a"],e=row["e"],
            inc=math.radians(row["i"]),
            Omega=math.radians(row["Omega"]),
            omega=math.radians(row["omega"]),
            M=math.radians(row["M"]),
            hash=row["designation"]
        )
    sim.move_to_com()
    sim.integrator="whfast"
    sim.dt=DT_YR
    sim.ri_whfast.safe_mode=0
    sim.ri_whfast.corrector=11
    return sim,start_idx

def run_one(tnos,a_p,e_p,Mdeg):
    m=mass_on_mu(a_p,e_p)
    sim,start=make_sim(tnos,a_p,e_p,m,Mdeg)
    initial=phase_order(sim,start)
    sim.integrate(T_YR,exact_finish_time=0)
    final=phase_order(sim,start)
    return {
        "phase_deg":Mdeg,"m_earth":m,
        "R_initial":initial["R"],"R_final":final["R"],
        "delta_R":final["R"]-initial["R"],
        "mean_varpi_initial_deg":initial["mean_varpi_deg"],
        "mean_varpi_final_deg":final["mean_varpi_deg"],
        "bound_final":final["bound"],
        "min_q_final_au":final["min_q_au"],
    }

def aggregate(rows):
    return {
        "mean_delta_R":sum(r["delta_R"] for r in rows)/len(rows),
        "min_delta_R":min(r["delta_R"] for r in rows),
        "max_delta_R":max(r["delta_R"] for r in rows),
        "mean_bound_final":sum(r["bound_final"] for r in rows)/len(rows),
        "min_q_final_au":min(r["min_q_final_au"] for r in rows if r["min_q_final_au"] is not None),
    }

def main():
    source_url,tnos=fetch_tnos()
    baseline=[]
    # no-perturber baseline represented by vanishing PBH mass at a dummy orbit
    for ph in PHASES_DEG:
        sim,start=make_sim(tnos,300.0,0.0,0.0,ph)
        initial=phase_order(sim,start)
        sim.integrate(T_YR,exact_finish_time=0)
        final=phase_order(sim,start)
        baseline.append({
          "phase_deg":ph,"R_initial":initial["R"],"R_final":final["R"],
          "delta_R":final["R"]-initial["R"],"bound_final":final["bound"],
          "min_q_final_au":final["min_q_au"]
        })
    candidates=[]
    for a in A_GRID:
        for e in E_GRID:
            runs=[run_one(tnos,a,e,ph) for ph in PHASES_DEG]
            agg=aggregate(runs)
            candidates.append({
              "a_au":a,"e":e,"m_earth":mass_on_mu(a,e),
              "runs":runs,**agg,
            })
            print("CANDIDATE",a,e,mass_on_mu(a,e),agg["mean_delta_R"],flush=True)
    base=aggregate(baseline)
    for c in candidates:
        c["delta_R_over_baseline"]=c["mean_delta_R"]-base["mean_delta_R"]
        c["score"]=c["delta_R_over_baseline"] - 0.05*max(0.0,20.0-c["mean_bound_final"])
    ranked=sorted(candidates,key=lambda x:x["score"],reverse=True)
    receipt={
      "schema":"PBH_REBOUND_RELATIONAL_STABILITY_PILOT_V0_1",
      "status":"ANALYZED",
      "source_url":source_url,
      "n_tnos":len(tnos),
      "integration":{"duration_yr":T_YR,"dt_yr":DT_YR,"integrator":"WHFast","phase_samples_deg":list(PHASES_DEG)},
      "fixed_relations":{"varpi_p_deg":VARPI_P_DEG,"mu_rel_mearth_per_au3":MU_REL},
      "baseline":base,
      "candidates":ranked,
      "best_candidate":ranked[0],
      "method_boundary":[
        "pilot tests short-term preservation/contraction of observed apsidal order, not formation from an unbiased primordial disk",
        "giant planets are circular/coplanar approximations",
        "candidate inclination fixed to zero in v0.1",
        "200 kyr is much shorter than published Gyr-scale stability experiments",
        "result is a ranking diagnostic, not a posterior or detection"
      ],
      "authority":"CANDIDATE_ONLY","canon_allowed":False
    }
    print("PBH_REBOUND_PILOT_SUMMARY="+json.dumps({
      "baseline":base,
      "ranking":[{k:c[k] for k in ("a_au","e","m_earth","mean_delta_R","delta_R_over_baseline","mean_bound_final","min_q_final_au","score")} for c in ranked]
    },sort_keys=True,separators=(",",":")),flush=True)
    print("PBH_REBOUND_PILOT_RECEIPT="+json.dumps(receipt,sort_keys=True,separators=(",",":")),flush=True)
    return 0

if __name__=="__main__":
    raise SystemExit(main())
