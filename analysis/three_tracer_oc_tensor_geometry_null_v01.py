#!/usr/bin/env python3
from __future__ import annotations

import json
import math
import numpy as np

SEED = 20260928
N = 200_000
ARCSEC_TO_RAD = math.pi/(180.0*3600.0)

PROFILES = {
    "all": {
        "VP": {"ra":51.740502024,"dec":2.264706711,"v":[-1.1786631749570087e-05,2.2615152683874284e-04]},
        "TG": {"ra":4.345899408,"dec":13.791573111,"v":[-2.6247696089711193e-05,-1.247966479691484e-05]},
        "Sed": {"ra":56.75455475,"dec":7.578368673,"v":[6.806180461773833e-05,-1.9622805099622816e-04]},
    },
    "gaia": {
        "VP": {"ra":52.541159515,"dec":2.485277679,"v":[-6.961334931244616e-04,1.8428030163705406e-04]},
        "TG": {"ra":3.177041596,"dec":12.978306056,"v":[-3.422434771802735e-05,-4.510116347653832e-05]},
        "Sed": {"ra":58.717750107,"dec":8.182263317,"v":[-1.9724662276312152e-04,-7.370083810209547e-04]},
    },
}
KEYS=["VP","TG","Sed"]

OBS = {
    "all_shape":0.01946748887590975,
    "all_split":0.06745014210166754,
    "gaia_shape":0.106563140980909,
    "gaia_split":0.371259523786135,
    "axis_overlap":0.9987214635162807,
    "tensor_cosine":0.9900345451483645,
    "all_abs_axis_dec":75.96287632643228,
    "gaia_abs_axis_dec":76.17119484994835,
}


def design(prof):
    H=[]
    for k in KEYS:
        d=prof[k]
        a=math.radians(d["ra"]); de=math.radians(d["dec"])
        n=np.array([math.cos(de)*math.cos(a),math.cos(de)*math.sin(a),math.sin(de)])
        era=np.array([-math.sin(a),math.cos(a),0.0])
        ed=np.array([-math.cos(a)*math.sin(de),-math.sin(a)*math.sin(de),math.cos(de)])
        for e in (era,ed):
            H.append([
                e[0]*n[0]-e[2]*n[2],
                e[0]*n[1]+e[1]*n[0],
                e[0]*n[2]+e[2]*n[0],
                e[1]*n[1]-e[2]*n[2],
                e[1]*n[2]+e[2]*n[1],
            ])
    return np.asarray(H,float)


def tensor_batch(t):
    out=np.empty((len(t),3,3),float)
    out[:,0,0]=t[:,0]
    out[:,0,1]=out[:,1,0]=t[:,1]
    out[:,0,2]=out[:,2,0]=t[:,2]
    out[:,1,1]=t[:,3]
    out[:,1,2]=out[:,2,1]=t[:,4]
    out[:,2,2]=-t[:,0]-t[:,3]
    return out


def best_metrics(T):
    vals,vecs=np.linalg.eigh(T)
    fro=np.sqrt(np.sum(vals*vals,axis=1))
    shapes=[];splits=[];mus=[]
    for k in range(3):
        o=[j for j in range(3) if j!=k]
        diff=np.abs(vals[:,o[1]]-vals[:,o[0]])
        shapes.append((diff/math.sqrt(2.0))/np.maximum(fro,1e-300))
        mu=vals[:,k]/2.0
        mus.append(mu)
        splits.append(diff/np.maximum(np.abs(mu),1e-300))
    shapes=np.stack(shapes,axis=1)
    splits=np.stack(splits,axis=1)
    mus=np.stack(mus,axis=1)
    idx=np.argmin(shapes,axis=1)
    rr=np.arange(len(T))
    return (
        shapes[rr,idx],
        splits[rr,idx],
        mus[rr,idx],
        vecs[rr,:,idx],
    )


def pct(mask):
    return float(np.mean(mask))


def main():
    rng=np.random.default_rng(SEED)
    La=np.linalg.pinv(design(PROFILES["all"]))
    Lg=np.linalg.pinv(design(PROFILES["gaia"]))
    va=np.asarray([PROFILES["all"][k]["v"] for k in KEYS],float)*ARCSEC_TO_RAD
    vg=np.asarray([PROFILES["gaia"][k]["v"] for k in KEYS],float)*ARCSEC_TO_RAD

    theta=rng.uniform(0.0,2.0*math.pi,size=(N,3))
    ct=np.cos(theta); st=np.sin(theta)
    Ya=np.empty((N,6));Yg=np.empty((N,6))
    for j in range(3):
        x,y=va[j]
        Ya[:,2*j]=ct[:,j]*x-st[:,j]*y
        Ya[:,2*j+1]=st[:,j]*x+ct[:,j]*y
        x,y=vg[j]
        Yg[:,2*j]=ct[:,j]*x-st[:,j]*y
        Yg[:,2*j+1]=st[:,j]*x+ct[:,j]*y

    Ta=tensor_batch(Ya@La.T)
    Tg=tensor_batch(Yg@Lg.T)
    sha,spa,mua,axa=best_metrics(Ta)
    shg,spg,mug,axg=best_metrics(Tg)

    overlap=np.abs(np.sum(axa*axg,axis=1))
    cos=np.sum(Ta*Tg,axis=(1,2))/(
        np.sqrt(np.sum(Ta*Ta,axis=(1,2)))*np.sqrt(np.sum(Tg*Tg,axis=(1,2)))
    )
    deca=np.degrees(np.arcsin(np.clip(axa[:,2],-1.0,1.0)))
    decg=np.degrees(np.arcsin(np.clip(axg[:,2],-1.0,1.0)))

    all_gate=(sha<=0.20)&(spa<=0.40)
    gaia_gate=(shg<=0.20)&(spg<=0.40)
    strict_shape=(sha<=OBS["all_shape"])&(shg<=OBS["gaia_shape"])
    strict=(
        strict_shape
        &(overlap>=OBS["axis_overlap"])
        &(cos>=OBS["tensor_cosine"])
    )
    receipt={
        "schema":"PBH_OC_TENSOR_COUPLED_GEOMETRY_NULL_V0_1",
        "seed":SEED,
        "draws":N,
        "null":"independent uniform tangent-plane rotation per object; same rotation applied to all-data and Gaia-family profile for that object; curvature magnitudes and cross-profile relative vector geometry preserved",
        "observed":OBS,
        "rates":{
            "all_sign_free_gate":pct(all_gate),
            "gaia_sign_free_gate":pct(gaia_gate),
            "both_sign_free_gate":pct(all_gate&gaia_gate),
            "all_shape_at_least_as_extreme":pct(sha<=OBS["all_shape"]),
            "gaia_shape_at_least_as_extreme":pct(shg<=OBS["gaia_shape"]),
            "both_shapes_at_least_as_extreme":pct(strict_shape),
            "axis_overlap_at_least_observed":pct(overlap>=OBS["axis_overlap"]),
            "tensor_cosine_at_least_observed":pct(cos>=OBS["tensor_cosine"]),
            "joint_strict_observed":pct(strict),
            "joint_strict_observed_and_both_negative_mu":pct(strict&(mua<0)&(mug<0)),
            "both_abs_axis_dec_at_least_observed":pct(
                (np.abs(deca)>=OBS["all_abs_axis_dec"])
                &(np.abs(decg)>=OBS["gaia_abs_axis_dec"])
            ),
        },
        "axis_geometry":{
            "all_abs_dec_quantiles_deg":{
                "q50":float(np.quantile(np.abs(deca),0.50)),
                "q90":float(np.quantile(np.abs(deca),0.90)),
                "q95":float(np.quantile(np.abs(deca),0.95)),
            },
            "gaia_abs_dec_quantiles_deg":{
                "q50":float(np.quantile(np.abs(decg),0.50)),
                "q90":float(np.quantile(np.abs(decg),0.90)),
                "q95":float(np.quantile(np.abs(decg),0.95)),
            },
        },
        "interpretation_boundary":{
            "post_hoc_statistic":True,
            "not_a_detection_p_value":True,
            "profiles_are_correlated":True,
            "jpl_fit_not_independent_of_mpc_astrometry":True,
            "geometry_null_only":True,
            "PBH_not_inferred":True,
        },
        "canon_allowed":False,
    }
    print("PBH_OC_TENSOR_NULL="+json.dumps(receipt,sort_keys=True),flush=True)


if __name__=="__main__":
    main()
