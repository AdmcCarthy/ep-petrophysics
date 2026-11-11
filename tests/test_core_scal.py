
import numpy as np
import pandas as pd

from ep_petrophysics import core_analysis as ca
from ep_petrophysics import core_scal as cs

LAB = {"low_max": 2.0, "low_a": 0.68, "low_b": 1.06, "m_a": 0.777, "m_b": 0.61, "iterations": 2}


def test_lab_klinkenberg_formula_low_and_high_range():
    kn2 = pd.Series([1.0, 13.8])
    kl = ca.lab_klinkenberg(kn2, pd.Series([0.5, 0.662]), LAB)
    assert abs(kl.iloc[0] - 0.68) < 1e-9                       # below 2 mD: Kl = 0.68 x Kn2^1.06
    assert abs(kl.iloc[1] - 11.5) < 0.2                        # a delivered plug of the final report: Kn2 13.8, 1/Pm 0.662 -> Kl 11.5
    assert np.isnan(ca.lab_klinkenberg(pd.Series([5.0]), pd.Series([np.nan]), LAB).iloc[0])


def test_klinkenberg_reproduction_flags_a_value_that_does_not_follow_the_formula():
    p = pd.DataFrame({"well_name": ["W", "W"], "depth_md_m": [1.0, 2.0], "kg_h_md": [13.8, 13.8], "inv_pm_h": [0.662, 0.662], "kl_h_md": [11.5, 5.0],
                      "kg_v_md": [np.nan, np.nan], "inv_pm_v": [np.nan, np.nan], "kl_v_md": [np.nan, np.nan]})
    r = ca.klinkenberg_reproduction(p, LAB, 0.02)
    assert r["agrees"].tolist() == [True, False]


def test_family_and_krw_check_flag():
    raw = pd.DataFrame({"depth_md_m": [3845.0, 3846.0], "kl_md": [100.0, 100.0], "ko_swi_md": [80.0, 80.0], "kw_sor_md": [80.0, 160.0], "krw_sor": [1.0, 1.0],
                        "wwi": [np.nan, 0.2], "owi": [np.nan, 0.1], "test_type": ["Cent., o/w rel.perm., imb. water", "USBM-test, restored"]})
    d, ev = cs.clean_endpoints(raw, "W", "f.xls", "s")
    assert d["test_family"].tolist() == ["centrifuge, imbibition oil relative permeability", "wettability, restored"]
    assert d["flag_krw_differs"].tolist() == [False, True]                       # 160 / 80 = 2.0 against the stated 1.0
    assert abs(d["wi_amott"].iloc[1] - 0.1) < 1e-9 and any(e["rule"] == "krw_sor_differs_from_kw_over_ko" for e in ev)


def test_waterflood_check_flags_a_copied_starting_point():
    ep = pd.DataFrame({"depth_md_m": [3911.0], "test_family": ["water flood, reservoir conditions"], "swi": [0.17], "sor": [0.24]})
    cv = pd.DataFrame({"nominal_depth_md_m": [3911.0] * 3, "sw": [0.088, 0.5, 0.76]})
    out, ev = cs.waterflood_checks(cv, ep, "s")
    assert out["flag_swi_mismatch"].iloc[0] and not out["flag_sor_mismatch"].iloc[0]


def test_corey_fit_recovers_the_exponent():
    sw = np.linspace(0.2, 0.8, 40)
    kro = ((1 - sw - 0.2) / (1 - 0.2 - 0.2)) ** 3.0
    curve = pd.DataFrame({"depth_md_m": 1.0, "sw": sw, "kro": kro})
    ep = pd.DataFrame({"depth_md_m": [1.0], "corey_oil": [4.0], "swi": [0.2], "sor": [0.2]})
    f = cs.fit_corey_oil(curve, ep)
    assert abs(f["corey_oil_fit"].iloc[0] - 3.0) < 0.05


