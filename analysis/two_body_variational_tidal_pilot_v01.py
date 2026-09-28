#!/usr/bin/env python3
from __future__ import annotations

import csv, hashlib, io, json, math
from collections import defaultdict
from datetime import datetime, timezone

import numpy as np
import requests

MPC_URL="https://data.minorplanetcenter.net/api/get-obs"
HORIZONS_URL="https://ssd.jpl.nasa.gov/api/horizons.api"
EXPECTED_ROWS_SHA256="78492d3031d0e3fd41c883a14117c79c66e5c165d072b2c0293749c12d9ebbf3"
ARCSEC_TO_RAD=math.pi/(180.0*3600.0)
DAY_PER_YEAR=365.25
MU_SUN_AU3_YR2=0.0002959122082855911*(DAY_PER_YEAR**2)
HMAX_YR=2.0/DAY_PER_YEAR
CHUNK=40

OBJECTS={
 "2012 VP113":"'DES=2012 VP113;'",
 "2015 TG387":"'541132;'",
 "Sedna":"'90377;'",
}

QBASIS=np.array([
 [[1,0,0],[0,0,0],[0,0,-1]],
 [[0,1,0],[1,0,0],[0,0,0]],
 [[0,0,1],[0,0,0],[1,0,0]],
 [[0,0,0],[0,1,0],[0,0,-1]],
 [[0,0,0],[0,0,1],[0,1,0]],
],dtype=float)


def f(x):
    try:
        v=float(x)
        return v if math.isfinite(v) else None
    except Exception:
        return None


def horizons_time(s):
    dt=datetime.fromisoformat(s.replace("Z","+00:00")).astimezone(timezone.utc)
    return dt.strftime("%Y-%b-%d %H:%M:%S.")+f"{dt.microsecond//1000:03d}"


def parse_observer(result):
    body=result.split("$$SOE",1)[1].split("$$EOE",1)[0].strip()
    rows=[]
    for rec in csv.reader(io.StringIO(body)):
        v=[x.strip() for x in rec]
        if len(v)<6: continue
        jd,ra,dec=f(v[1]),f(v[4]),f(v[5])
        if jd is None or ra is None or dec is None:
            raise RuntimeError("observer parse:"+repr(v))
        rows.append((jd,ra,dec))
    return rows


def parse_vectors(result):
    body=result.split("$$SOE",1)[1].split("$$EOE",1)[0].strip()
    rows=[]
    for rec in csv.reader(io.StringIO(body)):
        v=[x.strip() for x in rec]
        if len(v)<8: continue
        nums=[f(v[i]) for i in (0,2,3,4,5,6,7)]
        if any(x is None for x in nums):
            raise RuntimeError("vector parse:"+repr(v))
        rows.append({
          "jd":nums[0],
          "r":np.asarray(nums[1:4],dtype=float),
          "v_day":np.asarray(nums[4:7],dtype=float),
        })
    return rows


def fetch_mpc(desig):
    r=requests.get(MPC_URL,json={"desigs":[desig],"output_format":["ADES_DF"],"ades_version":"2022"},timeout=90)
    r.raise_for_status()
    p=r.json()
    if not p or "ADES_DF" not in p[0]: raise RuntimeError(desig+":ADES_DF missing")
    out=[]
    for row in p[0]["ADES_DF"]:
        ra,dec=f(row.get("ra")),f(row.get("dec"))
        stn=str(row.get("stn") or "").strip()
        tm=str(row.get("obstime") or "").strip()
        if ra is None or dec is None or not stn or not tm: continue
        out.append({
          "obs_id":row.get("obsid"),"time":tm,"stn":stn,
          "catalog":str(row.get("astcat") or "UNKNOWN"),
          "ra_obs_deg":ra,"dec_obs_deg":dec,
        })
    return out


def observer_chunk(command,stn,rows):
    tlist=" ".join("'"+horizons_time(x["time"])+"'" for x in rows)
    p={"format":"json","COMMAND":command,"OBJ_DATA":"'NO'","MAKE_EPHEM":"'YES'",
       "EPHEM_TYPE":"'OBSERVER'","CENTER":f"'{stn}@399'","TLIST":tlist,
       "TLIST_TYPE":"'CAL'","TIME_TYPE":"'UT'","QUANTITIES":"'1'",
       "REF_SYSTEM":"'ICRF'","ANG_FORMAT":"'DEG'","CSV_FORMAT":"'YES'",
       "EXTRA_PREC":"'YES'","CAL_FORMAT":"'BOTH'"}
    r=requests.get(HORIZONS_URL,params=p,timeout=90); r.raise_for_status()
    result=r.json().get("result","")
    if "$$SOE" not in result or "** Error" in result: raise RuntimeError(result[:700])
    z=parse_observer(result)
    if len(z)!=len(rows): raise RuntimeError(f"{stn}:observer count {len(z)} != {len(rows)}")
    return z


def fetch_oc(desig,command):
    obs=fetch_mpc(desig); by=defaultdict(list)
    for x in obs: by[x["stn"]].append(x)
    out=[]; failures=[]
    for stn in sorted(by):
        rr=sorted(by[stn],key=lambda x:x["time"])
        try:
            pred=[]
            for i in range(0,len(rr),CHUNK):
                pred.extend(observer_chunk(command,stn,rr[i:i+CHUNK]))
            for o,(jd,ra_p,dec_p) in zip(rr,pred,strict=True):
                dra=((o["ra_obs_deg"]-ra_p+180.0)%360.0)-180.0
                out.append({
                  "object":desig,"obs_id":o["obs_id"],"time":o["time"],"jd_ut":jd,
                  "stn":stn,"catalog":o["catalog"],
                  "ra_obs_deg":o["ra_obs_deg"],"dec_obs_deg":o["dec_obs_deg"],
                  "ra_pred_deg":ra_p,"dec_pred_deg":dec_p,
                  "dra_cosdec_arcsec":dra*3600.0*math.cos(math.radians(dec_p)),
                  "ddec_arcsec":(o["dec_obs_deg"]-dec_p)*3600.0,
                  "rho_arcsec":math.hypot(dra*3600.0*math.cos(math.radians(dec_p)),
                                          (o["dec_obs_deg"]-dec_p)*3600.0),
                })
        except Exception as exc:
            failures.append({"stn":stn,"n":len(rr),"error":type(exc).__name__+":"+str(exc)})
    out.sort(key=lambda x:(x["jd_ut"],x["stn"],str(x["obs_id"])))
    return out,failures


def nightly(rows):
    g=defaultdict(list)
    for r in rows:
        key=(r["stn"],r["catalog"],r["time"][:10])
        g[key].append(r)
    out=[]
    for (stn,cat,date),rr in sorted(g.items(),key=lambda kv:min(x["jd_ut"] for x in kv[1])):
        out.append({
          "stn":stn,"catalog":cat,"date":date,"n_raw":len(rr),
          "jd_ut":sum(x["jd_ut"] for x in rr)/len(rr),
          "dra_arcsec":sum(x["dra_cosdec_arcsec"] for x in rr)/len(rr),
          "ddec_arcsec":sum(x["ddec_arcsec"] for x in rr)/len(rr),
        })
    out.sort(key=lambda x:x["jd_ut"])
    return out


def vector_query(command,center,jds):
    out=[]
    for i in range(0,len(jds),CHUNK):
        chunk=jds[i:i+CHUNK]
        tlist=" ".join(f"'{x:.12f}'" for x in chunk)
        p={"format":"json","COMMAND":command,"OBJ_DATA":"'NO'","MAKE_EPHEM":"'YES'",
           "EPHEM_TYPE":"'VECTORS'","CENTER":center,"TLIST":tlist,"TLIST_TYPE":"'JD'",
           "TIME_TYPE":"'UT'","REF_SYSTEM":"'ICRF'","REF_PLANE":"'FRAME'",
           "OUT_UNITS":"'AU-D'","VEC_TABLE":"'2'","CSV_FORMAT":"'YES'"}
        r=requests.get(HORIZONS_URL,params=p,timeout=90); r.raise_for_status()
        result=r.json().get("result","")
        if "$$SOE" not in result or "** Error" in result: raise RuntimeError(result[:700])
        out.extend(parse_vectors(result))
    if len(out)!=len(jds): raise RuntimeError(f"vector count {len(out)} != {len(jds)}")
    return out


def grav_D(r):
    rr=float(r@r); R=math.sqrt(rr)
    I=np.eye(3)
    return -MU_SUN_AU3_YR2*(I/(R**3)-3.0*np.outer(r,r)/(R**5))


def acc(r):
    R=float(np.linalg.norm(r))
    return -MU_SUN_AU3_YR2*r/(R**3)


# Y layout:
# target r3,v3, Phi36, Stide30, earth r3,v3, Etide30
def unpack(Y):
    k=0
    rt=Y[k:k+3]; k+=3
    vt=Y[k:k+3]; k+=3
    Phi=Y[k:k+36].reshape(6,6); k+=36
    St=Y[k:k+30].reshape(6,5); k+=30
    re=Y[k:k+3]; k+=3
    ve=Y[k:k+3]; k+=3
    Se=Y[k:k+30].reshape(6,5); k+=30
    return rt,vt,Phi,St,re,ve,Se


def pack(rt,vt,Phi,St,re,ve,Se):
    return np.concatenate([rt,vt,Phi.ravel(),St.ravel(),re,ve,Se.ravel()])


def deriv(Y):
    rt,vt,Phi,St,re,ve,Se=unpack(Y)
    Dt=grav_D(rt); De=grav_D(re)
    Jt=np.block([[np.zeros((3,3)),np.eye(3)],[Dt,np.zeros((3,3))]])
    Je=np.block([[np.zeros((3,3)),np.eye(3)],[De,np.zeros((3,3))]])
    Ft=np.column_stack([Q@rt for Q in QBASIS])
    Fe=np.column_stack([Q@re for Q in QBASIS])
    dSt=Jt@St + np.vstack([np.zeros((3,5)),Ft])
    dSe=Je@Se + np.vstack([np.zeros((3,5)),Fe])
    return pack(vt,acc(rt),Jt@Phi,dSt,ve,acc(re),dSe)


def rk4(Y,h):
    k1=deriv(Y)
    k2=deriv(Y+0.5*h*k1)
    k3=deriv(Y+0.5*h*k2)
    k4=deriv(Y+h*k3)
    return Y+(h/6.0)*(k1+2*k2+2*k3+k4)


def integrate_offsets(Y0,offsets):
    result={}
    pos=sorted([x for x in offsets if x>=0])
    neg=sorted([x for x in offsets if x<0],reverse=True)
    for targets in (pos,neg):
        Y=Y0.copy(); t=0.0
        for target in targets:
            direction=1.0 if target>=t else -1.0
            while abs(target-t)>1e-13:
                h=direction*min(HMAX_YR,abs(target-t))
                Y=rk4(Y,h); t+=h
            result[target]=Y.copy()
    return result


def tangent_from_vec(r):
    rho=float(np.linalg.norm(r)); u=r/rho
    ra=math.atan2(u[1],u[0]); dec=math.asin(u[2])
    ea=np.array([-math.sin(ra),math.cos(ra),0.0])
    ed=np.array([-math.cos(ra)*math.sin(dec),-math.sin(ra)*math.sin(dec),math.cos(dec)])
    return rho,np.vstack([ea,ed])


def normalize_columns(X):
    norms=np.linalg.norm(X,axis=0)
    keep=norms>1e-14
    return X[:,keep]/norms[keep],keep,norms


def nuisance_matrix(nights,A):
    n=len(nights); blocks=[A]
    # coordinate intercepts
    cra=np.zeros(2*n); cdec=np.zeros(2*n)
    cra[0::2]=1.0; cdec[1::2]=1.0
    blocks.extend([cra[:,None],cdec[:,None]])
    for key in ("stn","catalog"):
        levels=sorted({x[key] for x in nights})
        for level in levels[1:]:
            ra=np.zeros(2*n); de=np.zeros(2*n)
            for i,x in enumerate(nights):
                if x[key]==level:
                    ra[2*i]=1.0; de[2*i+1]=1.0
            blocks.extend([ra[:,None],de[:,None]])
    return np.column_stack(blocks)


def object_design(name,command,raw,nights):
    jds=[x["jd_ut"] for x in nights]
    tref=sum(jds)/len(jds)
    # heliocentric target and Earth state at same reference epoch
    rt0=vector_query(command,"'500@10'",[tref])[0]
    re0=vector_query("'399'","'500@10'",[tref])[0]
    Y0=pack(
      rt0["r"],rt0["v_day"]*DAY_PER_YEAR,np.eye(6),np.zeros((6,5)),
      re0["r"],re0["v_day"]*DAY_PER_YEAR,np.zeros((6,5))
    )
    offsets=[(jd-tref)/DAY_PER_YEAR for jd in jds]
    states=integrate_offsets(Y0,offsets)
    geo=vector_query(command,"'500@399'",jds)

    A=np.zeros((2*len(nights),6)); B=np.zeros((2*len(nights),5)); y=np.zeros(2*len(nights))
    for i,(night,g) in enumerate(zip(nights,geo,strict=True)):
        Y=states[offsets[i]]
        rt,vt,Phi,St,re,ve,Se=unpack(Y)
        rho,E=tangent_from_vec(g["r"])
        M=E/rho
        A[2*i:2*i+2,:]=M@Phi[:3,:]
        B[2*i:2*i+2,:]=M@(St[:3,:]-Se[:3,:])
        y[2*i]=night["dra_arcsec"]*ARCSEC_TO_RAD
        y[2*i+1]=night["ddec_arcsec"]*ARCSEC_TO_RAD

    X=nuisance_matrix(nights,A)
    Xn,keep,norms=normalize_columns(X)
    beta=np.linalg.lstsq(Xn,y,rcond=1e-11)[0]
    gamma=np.linalg.lstsq(Xn,B,rcond=1e-11)[0]
    yp=y-Xn@beta
    Bp=B-Xn@gamma
    s=np.linalg.svd(Bp,compute_uv=False)
    colfrac=[]
    for k in range(5):
        den=float(B[:,k]@B[:,k]); num=float(Bp[:,k]@Bp[:,k])
        colfrac.append(num/den if den>0 else 0.0)
    return {
      "t_ref_jd_ut":tref,"n_raw":len(raw),"n_nights":len(nights),
      "A":A,"B":B,"y":y,"X":X,"yp":yp,"Bp":Bp,
      "diag":{
        "rank_nuisance":int(np.linalg.matrix_rank(Xn)),
        "n_nuisance_columns_kept":int(Xn.shape[1]),
        "Bp_rank":int(np.linalg.matrix_rank(Bp)),
        "Bp_singular_values_yr2":s.tolist(),
        "Bp_condition":float(s[0]/s[-1]) if len(s) and s[-1]>0 else float("inf"),
        "tide_column_information_fraction_after_nuisance":colfrac,
        "y_rms_arcsec":float(np.sqrt(np.mean(y*y))/ARCSEC_TO_RAD),
        "yp_rms_arcsec":float(np.sqrt(np.mean(yp*yp))/ARCSEC_TO_RAD),
        "model_boundary":"two-body Sun variational dynamics; geometric geocentric projection; no planetary-gradient/light-time derivative",
      }
    }


def vec5_to_T(t):
    a,b,c,d,e=t
    return np.array([[a,b,c],[b,d,e],[c,e,-a-d]],dtype=float)


def main():
    raw_all=[]; raw_by={}; failures={}
    for name,command in OBJECTS.items():
        print("PBH_VAR_PROGRESS="+json.dumps({"stage":"oc","object":name}),flush=True)
        rows,fail=fetch_oc(name,command)
        raw_by[name]=rows; raw_all.extend(rows); failures[name]=fail
    canonical=json.dumps(raw_all,sort_keys=True,separators=(",",":")).encode()
    raw_sha=hashlib.sha256(canonical).hexdigest()
    if raw_sha!=EXPECTED_ROWS_SHA256:
        raise RuntimeError("frozen raw O-C hash mismatch:"+raw_sha)
    if any(failures[x] for x in failures):
        raise RuntimeError("station failures:"+json.dumps(failures,sort_keys=True))

    designs={}
    Bstack=[]; ystack=[]
    nightly_hash_input={}
    for name,command in OBJECTS.items():
        ns=nightly(raw_by[name])
        nightly_hash_input[name]=ns
        print("PBH_VAR_PROGRESS="+json.dumps({"stage":"variational","object":name,"n_nights":len(ns)}),flush=True)
        d=object_design(name,command,raw_by[name],ns)
        designs[name]=d
        Bstack.append(d["Bp"]); ystack.append(d["yp"])

    B=np.vstack(Bstack); y=np.concatenate(ystack)
    sv=np.linalg.svd(B,compute_uv=False)
    rank=int(np.linalg.matrix_rank(B))
    t=np.linalg.lstsq(B,y,rcond=1e-11)[0]
    pred=B@t
    rss0=float(y@y); rss=float((y-pred)@(y-pred))
    T=vec5_to_T(t)
    eigvals,eigvecs=np.linalg.eigh(T)
    imax=int(np.argmax(eigvals))
    axis=eigvecs[:,imax]
    mu=float(eigvals[imax]/2.0)
    # positive-point-mass eigenshape diagnostic
    Q=3*np.outer(axis,axis)-np.eye(3)
    mu_ls=float(np.sum(T*Q)/6.0)
    shape=float(np.linalg.norm(T-mu_ls*Q,"fro")/max(np.linalg.norm(T,"fro"),1e-300))
    trans=float(abs(eigvals[1]-eigvals[0])/max(abs(mu_ls),1e-300))
    ra=(math.degrees(math.atan2(axis[1],axis[0]))%360.0)
    dec=math.degrees(math.asin(axis[2]/np.linalg.norm(axis)))
    nightly_sha=hashlib.sha256(json.dumps(nightly_hash_input,sort_keys=True,separators=(",",":")).encode()).hexdigest()

    compact={
      "schema":"PBH_TWO_BODY_VARIATIONAL_TIDAL_PILOT_V0_1",
      "status":"ANALYZED_DIAGNOSTIC_ONLY",
      "raw_rows_sha256":raw_sha,"nightly_rows_sha256":nightly_sha,
      "objects":{k:v["diag"]|{"n_raw":v["n_raw"],"n_nights":v["n_nights"],"t_ref_jd_ut":v["t_ref_jd_ut"]} for k,v in designs.items()},
      "common":{
        "n_rows":int(B.shape[0]),"rank":rank,"singular_values_yr2":sv.tolist(),
        "condition":float(sv[0]/sv[-1]) if sv[-1]>0 else float("inf"),
        "tensor_vec5_yr2":t.tolist(),"tensor_3x3_yr2":T.tolist(),
        "eigenvalues_yr2":eigvals.tolist(),
        "axis_largest_eigenvalue":{"ra_deg":ra,"dec_deg":dec,"vector":axis.tolist()},
        "mu_largest_over2_yr2":mu,"mu_frobenius_yr2":mu_ls,
        "point_mass_shape_residual":shape,
        "transverse_split_over_mu":trans,
        "nuisance_profiled_rss0":rss0,"rss_tide":rss,
        "explained_fraction":1.0-rss/rss0 if rss0>0 else 0.0,
        "nuisance_profiled_rms_arcsec":float(np.sqrt(rss0/len(y))/ARCSEC_TO_RAD),
        "post_tide_rms_arcsec":float(np.sqrt(rss/len(y))/ARCSEC_TO_RAD),
      },
      "interpretation_boundary":{
        "not_full_N_body":True,"not_independent_of_JPL_fit":True,
        "no_planetary_gradient":True,"no_light_time_derivative":True,
        "nightly_equal_weight":True,"PBH_not_inferred":True,
        "purpose":"test whether physical tide response survives local-state plus station/catalog nuisance projection before high-fidelity N-body build",
      },
      "canon_allowed":False,
    }
    print("PBH_VAR_SUMMARY="+json.dumps(compact,sort_keys=True),flush=True)

if __name__=="__main__":
    main()
