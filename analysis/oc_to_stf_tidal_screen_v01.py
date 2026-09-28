#!/usr/bin/env python3
from __future__ import annotations
import csv, io, json, math
from datetime import datetime, timedelta, timezone

import numpy as np
import requests

HORIZONS="https://ssd.jpl.nasa.gov/api/horizons.api"
ARCSEC_TO_RAD=math.pi/(180.0*3600.0)

# Frozen from PBH_MPC_HORIZONS_OC_RESIDUALS_V0_2.
DATA={
 "2012 VP113":{
   "command":"'DES=2012 VP113;'",
   "t0_jd_ut":2457126.7521158545,
   "dra2_arcsec_yr2":-1.1786631749570087e-05,
   "ddec2_arcsec_yr2":2.2615152683874284e-04,
 },
 "2015 TG387":{
   "command":"'541132;'",
   "t0_jd_ut":2458704.6562090847,
   "dra2_arcsec_yr2":-2.6247696089711193e-05,
   "ddec2_arcsec_yr2":-1.2479664796914840e-05,
 },
 "Sedna":{
   "command":"'90377;'",
   "t0_jd_ut":2458049.6644799025,
   "dra2_arcsec_yr2":6.806180461773833e-05,
   "ddec2_arcsec_yr2":-1.9622805099622816e-04,
 },
}

def jd_to_cal(jd):
    dt=datetime(1970,1,1,tzinfo=timezone.utc)+timedelta(days=float(jd)-2440587.5)
    return dt.strftime("%Y-%b-%d %H:%M:%S.")+f"{dt.microsecond//1000:03d}"

def query_radec(command,jd):
    params={
      "format":"json","COMMAND":command,"OBJ_DATA":"'NO'","MAKE_EPHEM":"'YES'",
      "EPHEM_TYPE":"'OBSERVER'","CENTER":"'500@399'",
      "TLIST":f"'{jd_to_cal(jd)}'","TLIST_TYPE":"'CAL'","TIME_TYPE":"'UT'",
      "QUANTITIES":"'1'","REF_SYSTEM":"'ICRF'","ANG_FORMAT":"'DEG'",
      "CSV_FORMAT":"'YES'","EXTRA_PREC":"'YES'","CAL_FORMAT":"'BOTH'",
    }
    r=requests.get(HORIZONS,params=params,timeout=60); r.raise_for_status()
    result=r.json().get("result","")
    if "$$SOE" not in result or "$$EOE" not in result:
        raise RuntimeError(result[:1000])
    body=result.split("$$SOE",1)[1].split("$$EOE",1)[0].strip()
    rec=next(csv.reader(io.StringIO(body)))
    vals=[x.strip() for x in rec]
    return float(vals[4]),float(vals[5])

def basis(alpha_deg,dec_deg):
    a=math.radians(alpha_deg); d=math.radians(dec_deg)
    u=np.array([math.cos(d)*math.cos(a), math.cos(d)*math.sin(a), math.sin(d)])
    ea=np.array([-math.sin(a), math.cos(a), 0.0])
    ed=np.array([-math.cos(a)*math.sin(d),-math.sin(a)*math.sin(d),math.cos(d)])
    return u,ea,ed

def design_row(e,u):
    ex,ey,ez=e; ux,uy,uz=u
    return np.array([
      ex*ux-ez*uz,
      ex*uy+ey*ux,
      ex*uz+ez*ux,
      ey*uy-ez*uz,
      ey*uz+ez*uy,
    ],dtype=float)

def vec5_to_T(t):
    a,b,c,d,e=t
    return np.array([[a,b,c],[b,d,e],[c,e,-a-d]],dtype=float)

def T_to_vec5(T):
    return np.array([T[0,0],T[0,1],T[0,2],T[1,1],T[1,2]],dtype=float)

def q_tensor(ra_deg,dec_deg):
    a=math.radians(ra_deg); d=math.radians(dec_deg)
    n=np.array([math.cos(d)*math.cos(a),math.cos(d)*math.sin(a),math.sin(d)])
    return 3.0*np.outer(n,n)-np.eye(3),n

def direct_point_fit(A,y):
    best=None
    # coarse sphere grid
    for dec in np.arange(-88.0,88.0001,2.0):
      for ra in np.arange(0.0,360.0,2.0):
        Q,n=q_tensor(float(ra),float(dec))
        q=A@T_to_vec5(Q)
        qq=float(q@q)
        if qq<=0: continue
        mu=float(q@y/qq)
        if mu<=0: continue
        res=y-mu*q
        rss=float(res@res)
        if best is None or rss<best["rss"]:
            best={"ra_deg":float(ra),"dec_deg":float(dec),"mu_yr2":mu,"rss":rss}
    if best is None: return None
    # local refinement
    c=best
    for step,span in [(0.5,3.0),(0.1,0.8),(0.02,0.18)]:
      base_ra,base_dec=c["ra_deg"],c["dec_deg"]
      for dec in np.arange(max(-89.99,base_dec-span),min(89.99,base_dec+span)+1e-12,step):
        for ra_off in np.arange(-span,span+1e-12,step):
          ra=(base_ra+float(ra_off))%360.0
          Q,n=q_tensor(ra,float(dec))
          q=A@T_to_vec5(Q); qq=float(q@q)
          if qq<=0: continue
          mu=float(q@y/qq)
          if mu<=0: continue
          res=y-mu*q; rss=float(res@res)
          if rss<c["rss"]:
            c={"ra_deg":ra,"dec_deg":float(dec),"mu_yr2":mu,"rss":rss}
    return c

rows=[]; y=[]; sky={}
for name,d in DATA.items():
    ra,dec=query_radec(d["command"],d["t0_jd_ut"])
    u,ea,ed=basis(ra,dec)
    rows.append(design_row(ea,u)); y.append(d["dra2_arcsec_yr2"]*ARCSEC_TO_RAD)
    rows.append(design_row(ed,u)); y.append(d["ddec2_arcsec_yr2"]*ARCSEC_TO_RAD)
    sky[name]={"t0_jd_ut":d["t0_jd_ut"],"ra_deg":ra,"dec_deg":dec,
               "curvature_arcsec_yr2":[d["dra2_arcsec_yr2"],d["ddec2_arcsec_yr2"]]}

A=np.vstack(rows); y=np.array(y)
t,resid,rank,svals=np.linalg.lstsq(A,y,rcond=None)
T=vec5_to_T(t)
pred=A@t; r=y-pred
rss0=float(y@y); rss=float(r@r)

vals,vecs=np.linalg.eigh(T)
axis=vecs[:,-1]
# deterministic sign only for display
for z in axis:
    if abs(z)>1e-15:
        if z<0: axis=-axis
        break
Q=3*np.outer(axis,axis)-np.eye(3)
mu=float(np.sum(T*Q)/6.0)
Tpoint=mu*Q
shape_res=float(np.linalg.norm(T-Tpoint,"fro")/max(np.linalg.norm(T,"fro"),1e-300))
trans_split=float(abs(vals[1]-vals[0])/max(abs(mu),1e-300))
axis_ra=(math.degrees(math.atan2(axis[1],axis[0]))%360.0)
axis_dec=math.degrees(math.asin(axis[2]/np.linalg.norm(axis)))

point=direct_point_fit(A,y)
if point:
    point["explained_fraction_vs_zero"]=1.0-point["rss"]/rss0
    point["rms_residual_arcsec_yr2"]=math.sqrt(point["rss"]/len(y))/ARCSEC_TO_RAD

earth_per_sun=332946.0487
screen_mass_over_r3_earth_per_au3=(mu/(4*math.pi*math.pi))*earth_per_sun

out={
 "schema":"PBH_OC_TO_STF_TIDAL_SCREEN_V0_1",
 "status":"DIAGNOSTIC_ONLY",
 "authority":"NONINDEPENDENT_OC_SCREEN",
 "source_receipt_schema":"PBH_MPC_HORIZONS_OC_RESIDUALS_V0_2",
 "source_rows_sha256":"78492d3031d0e3fd41c883a14117c79c66e5c165d072b2c0293749c12d9ebbf3",
 "sky_at_profile_epoch":sky,
 "linear_system":{
   "n_equations":int(A.shape[0]),"n_stf_parameters":int(A.shape[1]),"rank":int(rank),
   "singular_values":svals.tolist(),"condition_number":float(svals[0]/svals[-1]),
 },
 "stf_fit":{
   "tensor_vec5_yr2":t.tolist(),"tensor_3x3_yr2":T.tolist(),
   "eigenvalues_yr2":vals.tolist(),"axis_sign_degenerate":{"ra_deg":axis_ra,"dec_deg":axis_dec,"vector":axis.tolist()},
   "mu_frobenius_projection_yr2":mu,
   "formal_mass_over_r3_earth_per_au3_DIAGNOSTIC_ONLY":screen_mass_over_r3_earth_per_au3,
   "point_mass_shape_residual":shape_res,
   "transverse_eigenvalue_split_over_mu":trans_split,
   "rss_rad2_yr4":rss,"zero_rss_rad2_yr4":rss0,
   "explained_fraction_vs_zero":1.0-rss/rss0,
   "rms_residual_arcsec_yr2":math.sqrt(rss/len(y))/ARCSEC_TO_RAD,
 },
 "direct_positive_point_mass_fit":point,
 "interpretation_boundary":{
   "current_Horizons_orbits_are_not_independent_of_MPC_astrometry":True,
   "O_C_polynomial_curvature_is_not_a_variational_tidal_measurement":True,
   "mu_and_mass_distance_ratio_are_screening_values_only":True,
   "no_significance_from_six_equations":True,
   "PBH_not_inferred":True,
 },
 "canon_allowed":False,
}
print("PBH_STF_SCREEN="+json.dumps(out,sort_keys=True),flush=True)
