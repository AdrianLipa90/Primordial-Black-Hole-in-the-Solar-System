#!/usr/bin/env python3
from __future__ import annotations

import json, math, urllib.parse
import numpy as np
import requests

API="https://ssd-api.jpl.nasa.gov/sbdb_query.api"
CONSTRAINT={"AND":["a|GE|150","q|GE|40","e|LT|1","data_arc|GE|365","n_obs_used|GE|5"]}
FIELDS=["spkid","pdes","full_name","epoch","e","a","q","i","om","w","ma",
        "sigma_i","sigma_om","sigma_w","data_arc","n_obs_used","condition_code","source"]
QUALITY={"condition_code_max":3,"data_arc_min_days":2000,"n_obs_used_min":15}
SEED=20260928
N_NULL=20000
N_BOOT=20000

def f(x):
    return None if x in (None,"") else float(x)

def peri_vec(i_deg,om_deg,w_deg):
    i,om,w=map(math.radians,(i_deg,om_deg,w_deg))
    return np.array([
        math.cos(om)*math.cos(w)-math.sin(om)*math.sin(w)*math.cos(i),
        math.sin(om)*math.cos(w)+math.cos(om)*math.sin(w)*math.cos(i),
        math.sin(w)*math.sin(i)
    ],dtype=float)

def pole_vec(i_deg,om_deg):
    i,om=map(math.radians,(i_deg,om_deg))
    return np.array([
        math.sin(i)*math.sin(om),
        -math.sin(i)*math.cos(om),
        math.cos(i)
    ],dtype=float)

def norm(v):
    n=float(np.linalg.norm(v))
    return v/n if n>0 else v

def lonlat(v):
    v=norm(v)
    return (math.degrees(math.atan2(v[1],v[0]))%360.0,
            math.degrees(math.asin(max(-1,min(1,v[2])))))

def ecl_to_eq(v):
    eps=math.radians(23.439291111)
    x,y,z=v
    return np.array([x,y*math.cos(eps)-z*math.sin(eps),y*math.sin(eps)+z*math.cos(eps)])

def radec(v):
    q=norm(ecl_to_eq(v))
    return (math.degrees(math.atan2(q[1],q[0]))%360.0,
            math.degrees(math.asin(max(-1,min(1,q[2])))))

def angsep(a,b):
    a=norm(a); b=norm(b)
    return math.degrees(math.acos(max(-1,min(1,float(np.dot(a,b))))))

def circ_mean_deg(vals):
    a=np.radians(vals)
    z=np.mean(np.exp(1j*a))
    return float(np.degrees(np.angle(z))%360.0), float(abs(z))

params={"fields":",".join(FIELDS),"sb-kind":"a","sb-class":"TNO",
        "sb-cdata":json.dumps(CONSTRAINT,separators=(",",":")),
        "full-prec":"true","sort":"a"}
url=API+"?"+urllib.parse.urlencode(params)
r=requests.get(url,timeout=120,headers={"User-Agent":"GREMLIN-PBH-RELATIONAL-ORBIT/0.1"})
r.raise_for_status()
raw=r.json()
fields=raw["fields"]
rows=[dict(zip(fields,row)) for row in raw["data"]]

sel=[]
for row in rows:
    try: cc=int(row["condition_code"])
    except: continue
    if cc>QUALITY["condition_code_max"]: continue
    if float(row["data_arc"])<QUALITY["data_arc_min_days"]: continue
    if int(row["n_obs_used"])<QUALITY["n_obs_used_min"]: continue
    sel.append(row)

P=np.stack([peri_vec(f(x["i"]),f(x["om"]),f(x["w"])) for x in sel])
H=np.stack([pole_vec(f(x["i"]),f(x["om"])) for x in sel])
meanP=np.mean(P,axis=0); meanH=np.mean(H,axis=0)
Rperi=float(np.linalg.norm(meanP)); Rpole=float(np.linalg.norm(meanH))
axisP=norm(meanP); axisH=norm(meanH)

Q=np.mean(np.einsum("ni,nj->nij",P,P),axis=0)-np.eye(3)/3.0
evals,evecs=np.linalg.eigh(Q)
order=np.argsort(evals)[::-1]
evals=evals[order]; evecs=evecs[:,order]
quad_axis=norm(evecs[:,0])
if np.dot(quad_axis,axisP)<0: quad_axis=-quad_axis

varpi=[(f(x["om"])+f(x["w"]))%360.0 for x in sel]
varpi_mean,varpi_R=circ_mean_deg(varpi)

rng=np.random.default_rng(SEED)
null_R=np.empty(N_NULL)
for k in range(N_NULL):
    Ps=[]
    for x in sel:
        i=f(x["i"])
        om=float(rng.uniform(0,360)); w=float(rng.uniform(0,360))
        Ps.append(peri_vec(i,om,w))
    null_R[k]=np.linalg.norm(np.mean(Ps,axis=0))
p_add1=float((1+np.sum(null_R>=Rperi))/(N_NULL+1))

boot_dirs=[]
for k in range(N_BOOT):
    ix=rng.integers(0,len(sel),len(sel))
    v=norm(np.mean(P[ix],axis=0))
    if np.dot(v,axisP)<0: v=-v
    boot_dirs.append(v)
boot_dirs=np.stack(boot_dirs)
boot_sep=np.array([angsep(axisP,v) for v in boot_dirs])

old_ra=107.576646; old_dec=-25.907579
ra=math.radians(old_ra); dec=math.radians(old_dec)
old_eq=np.array([math.cos(dec)*math.cos(ra),math.cos(dec)*math.sin(ra),math.sin(dec)])
eps=math.radians(23.439291111)
x,y,z=old_eq
old_ecl=np.array([x,y*math.cos(eps)+z*math.sin(eps),-y*math.sin(eps)+z*math.cos(eps)])

details=[]
for x,p,h,vp in zip(sel,P,H,varpi):
    details.append({
      "pdes":x["pdes"],"full_name":x["full_name"].strip(),
      "a_au":f(x["a"]),"q_au":f(x["q"]),"e":f(x["e"]),
      "i_deg":f(x["i"]),"node_deg":f(x["om"]),"argperi_deg":f(x["w"]),
      "varpi_deg":vp,"data_arc_days":float(x["data_arc"]),
      "n_obs_used":int(x["n_obs_used"]),"condition_code":int(x["condition_code"]),
      "perihelion_unit_ecliptic":p.tolist(),"pole_unit_ecliptic":h.tolist()
    })

out={
 "schema":"PBH_RELATIONAL_ORBIT_GEOMETRY_V0_1",
 "status":"ANALYZED",
 "source":"JPL_SBDB_QUERY_API",
 "source_url":url,
 "selection_parent_count":len(rows),
 "quality_cut":QUALITY,
 "n_selected":len(sel),
 "selected":details,
 "apsidal":{
   "mean_resultant_length_3d":Rperi,
   "mean_axis_ecliptic_lonlat_deg":lonlat(axisP),
   "mean_axis_equatorial_radec_deg":radec(axisP),
   "quadrupole_eigenvalues":evals.tolist(),
   "quadrupole_axis_ecliptic_lonlat_deg":lonlat(quad_axis),
   "quadrupole_axis_equatorial_radec_deg":radec(quad_axis),
   "varpi_circular_mean_deg":varpi_mean,
   "varpi_resultant_length":varpi_R,
   "isotropic_orientation_null_add_one_rank":p_add1,
   "bootstrap_axis_sep_deg":{"median":float(np.median(boot_sep)),"p68":float(np.quantile(boot_sep,.68)),"p95":float(np.quantile(boot_sep,.95))}
 },
 "poles":{
   "mean_resultant_length_3d":Rpole,
   "mean_axis_ecliptic_lonlat_deg":lonlat(axisH),
   "mean_axis_equatorial_radec_deg":radec(axisH)
 },
 "old_2025_candidate_comparison":{
   "old_candidate_ra_dec_deg":[old_ra,old_dec],
   "angular_separation_to_mean_apsidal_axis_deg":angsep(old_ecl,axisP),
   "angular_separation_to_antipode_mean_apsidal_axis_deg":angsep(old_ecl,-axisP),
   "angular_separation_to_quadrupole_axis_deg":min(angsep(old_ecl,quad_axis),angsep(old_ecl,-quad_axis))
 },
 "interpretation_boundary":"ORIENTATION_GEOMETRY_ONLY_NOT_ASTROMETRIC_O_MINUS_C_NOT_PBH_POSITION",
 "authority":"CANDIDATE_ONLY",
 "canon_allowed":False
}
print("PBH_REL_ORBIT_SUMMARY="+json.dumps({
 "n_selected":out["n_selected"],
 "apsidal":out["apsidal"],
 "poles":out["poles"],
 "old_comparison":out["old_2025_candidate_comparison"],
 "selected_names":[x["pdes"] for x in details]
},sort_keys=True,separators=(",",":")),flush=True)
print("PBH_REL_ORBIT_RECEIPT="+json.dumps(out,sort_keys=True,separators=(",",":")),flush=True)
