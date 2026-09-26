#!/usr/bin/env python3
from __future__ import annotations
import json, math

G=6.67430e-11
C=299792458.0
AU=149597870700.0
PC=3.085677581491367e16
M_SUN=1.98847e30
M_EARTH=5.9722e24
D_L_AU=583.1377550476387
D_L=D_L_AU*AU
D_S=1000.0*PC
V_EARTH=30000.0
V_ORBIT=989.7283900996955

SCENARIOS={
    "origins_5_earth_mass":5.0*M_EARTH,
    "narrative_10_earth_mass":10.0*M_EARTH,
    "published_5e_5_Msun":5.0e-5*M_SUN,
}

def lensing(M):
    DLS=D_S-D_L
    theta=math.sqrt((4.0*G*M/C**2)*(DLS/(D_L*D_S)))
    theta_mas=theta*206265000.0
    re=theta*D_L
    rs=2.0*G*M/C**2
    centroid_max_mas=theta_mas/(2.0*math.sqrt(2.0))
    return {
        "mass_kg":M,
        "mass_Earth":M/M_EARTH,
        "schwarzschild_radius_m":rs,
        "schwarzschild_diameter_m":2.0*rs,
        "theta_E_mas_at_1kpc_source":theta_mas,
        "max_astrometric_centroid_shift_mas":centroid_max_mas,
        "physical_Einstein_radius_km":re/1000.0,
        "tE_minutes_at_30km_s":re/V_EARTH/60.0,
        "tE_hours_at_0p9897km_s":re/V_ORBIT/3600.0,
    }

def main():
    parallax_rad=math.asin(1.0/D_L_AU)
    parallax_arcmin=math.degrees(parallax_rad)*60.0
    annual_parallax_speed_arcsec_day=(
        2.0*math.pi*parallax_rad*206265.0/365.25
    )
    orbital_transverse_AU_yr=0.13981923831769566
    orbital_pm_arcsec_yr=orbital_transverse_AU_yr/D_L_AU*206265.0

    receipt={
      "schema":"PBH_MICROLENSING_OBSERVATION_GATE_V0_1",
      "status":"ANALYZED",
      "geometry":{
        "lens_distance_AU":D_L_AU,
        "assumed_background_source_distance_pc":1000.0,
        "annual_parallax_semi_amplitude_arcmin":parallax_arcmin,
        "approx_max_parallax_track_speed_arcsec_day":annual_parallax_speed_arcsec_day,
        "kepler_branch_transverse_proper_motion_arcsec_yr":orbital_pm_arcsec_yr,
      },
      "mass_scenarios":{k:lensing(v) for k,v in SCENARIOS.items()},
      "key_result":{
        "old_1_to_3_day_sampling":"INADEQUATE_FOR_PHOTOMETRIC_MICROLENSING_AT_THIS_DISTANCE",
        "reason":"Earth-induced relative motion can drive Einstein-radius crossing times of order 1.5-2.8 minutes for 5-16.65 Earth masses at about 583 AU.",
        "astrometric_channel":"More forgiving geometrically because centroid shifts are of order 2.3-4.2 mas at u~sqrt(2), but still needs a known source track and cadence matched to the event.",
      },
      "formulae":{
        "theta_E":"sqrt[(4GM/c^2) * D_ls/(D_l D_s)]",
        "R_E":"D_l * theta_E",
        "t_E":"R_E / v_perp",
        "delta_theta_max":"theta_E/(2 sqrt(2)) for a point lens",
        "annual_parallax":"asin(1 AU / D_l)",
      },
      "epistemic_boundary":"These are point-lens predictions conditional on the candidate distance and mass scenarios. They do not establish that the lens exists or that it is a PBH.",
      "authority":"OBSERVATION_DESIGN",
      "canon_allowed":False
    }
    print(json.dumps(receipt,indent=2,sort_keys=True))

if __name__=="__main__":
    main()
