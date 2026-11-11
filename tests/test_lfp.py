import numpy as np
import pandas as pd

from ep_petrophysics import lfp as pe

LAS = """~Version Information
 VERS. 2.0 : x
~Curve Information
#MNEM.UNIT API CODE Curve
 DEPTH .M : d
 LFP_AI .kPa.s/m : a
 LFP_RHOB .g/cm3 : r
 LFP_VP .m/s : v
~Parameter
~A DEPTH LFP_AI LFP_RHOB LFP_VP
 100.0 9000.0 2.5
 3600.0
 100.5 -999.25 2.4 3500.0
 -9999.0 0.0 0.0 0.0
"""

DOC = """** 15/9-19A *******
DATE: 03/2002
Interpreter: Test Person
UNITS:
      3821.50 - 3919.60 Hugin Fm.
      3919.60 - 3966.33 Sleipner Fm.
EVALUATION: 3655 - 4097 m
WATER DEPTH:  85 m (completion log)
 SHALE VOLUME:
Hugin Fm. |  LFP_GRMIN  |  LFP_GRMAX  | LFP_NPHISH | 2.00  | GR/DN
SALINITY:  150000 ppm
  "LFP_PRESS" = -32.88298797607 + 0.09235171974 * DEPTH
      API     25.89   (0.899)  [API] (g/cc)
      GOR       109      [Sm3/Sm3]
  gas gravity   0.84  (oil associated gas)
MATRIX DENSITY: set to 2.65 g/cc general
"""






CFG = {"checks": {"impedance_rel": 0.02, "vp_from_dt_rel": 0.02, "plausible": {"ai_min": 1000.0, "phie_max": 0.5, "gr_max": 300.0}}}


def _curve():
    md = np.arange(100.0, 160.0, 1.0)
    rho = np.full(len(md), 2.4)
    vp = np.full(len(md), 3500.0)
    c = pd.DataFrame({"well_name_file": "W", "md_m": md, "lfp_rhob": rho, "lfp_vp": vp, "lfp_ai": rho * vp, "lfp_dt": 304800.0 / vp, "lfp_vs": 2000.0, "lfp_vpvs": 3500.0 / 2000.0,
                      "lfp_baddata": 0.0, "lfp_phie": 0.2, "lfp_gr": 50.0, "lfp_sand": 1.0})
    c["lfp_vp_b"], c["lfp_vp_o"], c["lfp_vp_g"] = vp * 1.0, vp * 0.97, vp * 0.96
    c["lfp_ai_b"], c["lfp_ai_o"], c["lfp_ai_g"] = c["lfp_ai"], c["lfp_ai"] * 0.94, c["lfp_ai"] * 0.90
    for p, v in (("rhob", 2.4), ("vpvs", 1.75)):
        c[f"lfp_{p}_b"], c[f"lfp_{p}_o"], c[f"lfp_{p}_g"] = v, v * 0.97, v * 0.94
    return c


def test_physical_checks_pass_on_consistent_logs_and_count_implausible_values():
    c = _curve()
    c.loc[3, "lfp_gr"] = 900.0
    c.loc[4, "lfp_ai"] = 5000.0                                                    # breaks impedance = density x velocity once
    r = pe.physical_checks(c, CFG).iloc[0]
    assert r["gr_above_limit"] == 1 and abs(r["impedance_within_tol"] - 59 / 60) < 1e-9 and r["vp_from_dt_within_tol"] == 1.0 and r["vpvs_within_tol"] == 1.0


def test_fluid_ordering_and_effect_use_the_net_sand_of_the_interval_only():
    c = _curve()
    c.loc[c["md_m"] >= 140, "lfp_sand"] = 0.0                                      # shale below 140 m is not net
    iv = {"W": (100.0, 150.0)}
    o = pe.fluid_ordering(c, iv, min_n=5).iloc[0]
    assert o["net_sand_samples"] == 40 and o["vp_brine_ge_oil"] == 1.0 and o["ai_oil_ge_gas"] == 1.0
    e = pe.fluid_effect(c, iv, min_n=5)
    ai = e[(e["property"] == "P-impedance") & (e["case"] == "gas")]["median_change_vs_brine_pct"].iloc[0]
    assert abs(ai - (-10.0)) < 1e-6                                                # gas 10 % below the brine case
    assert pe.fluid_effect(c, iv, min_n=500)["median_change_vs_brine_pct"].isna().all()    # fewer samples than min_n: blank, never zero


def test_report_checks_compare_statements_with_curves():
    md = np.arange(3500.0, 3600.0, 1.0)
    phie = np.where(md > 3550, 0.2, 0.05)
    c = pd.DataFrame({"well_name_file": "W", "md_m": md, "lfp_phie": phie, "lfp_phit": phie + 0.06, "lfp_vsh": 0.1, "lfp_sand": (phie > 0.1).astype(float), "lfp_subst": (phie > 0.1).astype(float),
                      "lfp_model": np.where((md >= 3560) & (md <= 3580), 3.0, 1.0), "lfp_temp": -2.68 + 0.02779 * md, "lfp_press": -32.88 + 0.09235 * md})
    f = pd.DataFrame([{"well_name_file": "R", "evaluation_start_md_m": 3500.0, "evaluation_end_md_m": 3600.0, "sand_phi_limit": 0.1, "sand_vsh_limit": 0.5, "rkb_msl_m": 78.0, "temperature_intercept": -2.68, "temperature_slope_per_m": 0.02779,
                       "pressure_intercept": -32.88, "pressure_slope_per_m": 0.09235}])
    p = pd.DataFrame({"well_name_file": ["R"], "parameter": ["LFP_OIL"], "value_text": ["1:3560-3580 m MD | 0: otherwise"]})
    cfg = {"report_to_curve": {"R": "W"}, "report_check": {"interval_m": 5.0, "rule_share": 0.99, "temperature_degC": 1.0, "pressure_bar": 2.0}}
    r = pe.report_checks(c, f, p, cfg).set_index("check")
    assert r.loc["oil_interval_vs_model_flag", "status"] == "agrees" and r.loc["sand_flag_rule_phie", "status"] == "agrees" and r.loc["sand_flag_rule_phit", "status"] == "differs"
    assert "MD RKB" in r.loc["temperature_line", "unit"] and r.loc["temperature_line", "agreement"] < 1e-6 and r.loc["subst_equals_sand", "agreement"] == 1.0