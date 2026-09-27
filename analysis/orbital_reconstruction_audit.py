#!/usr/bin/env python3
from __future__ import annotations
import json, math

POSITION_AU=(-158.399,500.044,-254.785)
REPORTED_DISTANCE_AU=576.0
MASS_MSUN=5.0e-5
A_AU=430.0
E=0.78
I_DEG=27.0
INTEGRATION_YR=100000.0
AU_M=149597870700.0
YEAR_S=365.25*86400.0
GM_SUN_SI=1.32712440018e20
MU=GM_SUN_SI*YEAR_S**2/AU_M**3*(1.0+MASS_MSUN)
EARTH_MASSES_PER_SOLAR=332946.0487
TWOPI=2.0*math.pi

def norm(v): return math.sqrt(sum(x*x for x in v))
def dot(a,b): return sum(x*y for x,y in zip(a,b))
def cross(a,b): return (a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0])
def matmul(A,B): return tuple(tuple(sum(A[i][k]*B[k][j] for k in range(3)) for j in range(3)) for i in range(3))
def matvec(A,v): return tuple(sum(A[i][k]*v[k] for k in range(3)) for i in range(3))
def R1(t): return ((1.,0.,0.),(0.,math.cos(t),-math.sin(t)),(0.,math.sin(t),math.cos(t)))
def R3(t): return ((math.cos(t),-math.sin(t),0.),(math.sin(t),math.cos(t),0.),(0.,0.,1.))

def orbital_from_state(rvec,vvec):
    r=norm(rvec); v=norm(vvec); h=cross(rvec,vvec); hmag=norm(h)
    vxh=cross(vvec,h); evec=tuple(vxh[j]/MU-rvec[j]/r for j in range(3))
    ecc=norm(evec); energy=v*v/2.0-MU/r
    aa=-MU/(2.0*energy) if abs(energy)>1e-15 else math.inf
    inc=math.degrees(math.acos(max(-1.,min(1.,h[2]/hmag)))) if hmag else None
    return {"a_AU":aa,"e":ecc,"i_deg":inc,"specific_energy_AU2_yr2":energy,"speed_AU_yr":v}

def reconstruct_branches():
    x,y,z=POSITION_AU; r=norm(POSITION_AU); inc=math.radians(I_DEG)
    p=A_AU*(1.0-E*E); cosf=(p/r-1.0)/E
    if not -1.0<=cosf<=1.0: raise RuntimeError("position incompatible with a,e")
    f0=math.acos(cosf); sin_u=z/(r*math.sin(inc))
    if not -1.0<=sin_u<=1.0: raise RuntimeError("position incompatible with inclination")
    asin_u=math.asin(sin_u)
    u_values=(asin_u%TWOPI,(math.pi-asin_u)%TWOPI)
    lam=math.atan2(y,x)%TWOPI; branches=[]
    for u in u_values:
        alpha=math.atan2(math.sin(u)*math.cos(inc),math.cos(u))
        Omega=(lam-alpha)%TWOPI
        for f in (f0,(-f0)%TWOPI):
            omega=(u-f)%TWOPI
            rot=matmul(matmul(R3(Omega),R1(inc)),R3(omega))
            r_pf=(r*math.cos(f),r*math.sin(f),0.0)
            vscale=math.sqrt(MU/p)
            v_pf=(-vscale*math.sin(f),vscale*(E+math.cos(f)),0.0)
            rr=matvec(rot,r_pf); vv=matvec(rot,v_pf)
            branches.append({
                "Omega_deg":math.degrees(Omega),"omega_deg":math.degrees(omega),
                "true_anomaly_deg":math.degrees(f),"argument_of_latitude_deg":math.degrees(u),
                "position_AU":list(rr),"velocity_AU_yr":list(vv),
                "speed_AU_yr":norm(vv),"speed_km_s":norm(vv)*AU_M/YEAR_S/1000.0,
                "radial_velocity_AU_yr":dot(vv,rr)/r,
                "position_closure_error_AU":norm(tuple(rr[j]-POSITION_AU[j] for j in range(3))),
            })
    return branches

def main():
    x,y,z=POSITION_AU; r=norm(POSITION_AU)
    ra=math.degrees(math.atan2(y,x))%360.; dec=math.degrees(math.asin(z/r))
    p=A_AU*(1-E*E); q=A_AU*(1-E); Q=A_AU*(1+E)
    period=TWOPI*math.sqrt(A_AU**3/MU)
    speed=math.sqrt(MU*(2.0/r-1.0/A_AU))
    branches=reconstruct_branches()
    closest=min(branches,key=lambda b:abs(b["Omega_deg"]-180.)+abs(b["omega_deg"]-130.))
    receipt={
      "schema":"PBH_ORBITAL_RECONSTRUCTION_AUDIT_V0_1","status":"ANALYZED",
      "source_model":{"position_AU":list(POSITION_AU),"reported_distance_AU":REPORTED_DISTANCE_AU,"mass_Msun":MASS_MSUN,"a_AU":A_AU,"e":E,"i_deg":I_DEG,"reported_integration_years":INTEGRATION_YR},
      "direct_checks":{"cartesian_norm_AU":r,"distance_report_error_AU":r-REPORTED_DISTANCE_AU,"distance_report_error_fraction":r/REPORTED_DISTANCE_AU-1.0,"ra_deg_from_cartesian":ra,"dec_deg_from_cartesian":dec,"mass_Earth_from_5e_5_Msun":MASS_MSUN*EARTH_MASSES_PER_SOLAR,"perihelion_AU":q,"aphelion_AU":Q,"semilatus_rectum_AU":p,"two_body_period_yr":period,"cycles_in_100kyr":INTEGRATION_YR/period,"vis_viva_speed_AU_yr":speed,"vis_viva_speed_km_s":speed*AU_M/YEAR_S/1000.0,"max_annual_parallax_arcmin":math.degrees(math.asin(1.0/r))*60.0},
      "historical_code_variant_checks":{"state_x576_vy0p001_AUyr":orbital_from_state((576.,0.,0.),(0.,0.001,0.)),"state_x576_vy1_AUyr":orbital_from_state((576.,0.,0.),(0.,1.,0.))},
      "kepler_consistent_state_branches":branches,
      "geometric_clue":{"branch_near_round_angles":closest,"inference_only":"One exact-consistency branch lands extremely close to Omega=180 deg and omega=130 deg. These angles are not present in the recovered source, so this is a reconstruction clue only."},
      "authority":"CANDIDATE_ONLY","canon_allowed":False
    }
    print(json.dumps(receipt,indent=2,sort_keys=True))

if __name__=="__main__": main()
