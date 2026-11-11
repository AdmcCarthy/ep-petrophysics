import numpy as np
import pandas as pd

from ep_petrophysics import petrophysics as pe

REPORT = """PETROPHYSICS VOLVE 15/9-F-12
1.4 Permeability
Klogh = 10**(8* Phif-9*Vsh+2)   Hugin
Klogh = 10**(32*Phif-2*Vsh-3)   Sleipner
Phif = Phid + A*(Nphi-Phid) + B
Evaluation Parameters  15/9-F12
Parameter Draupne Heather Hugin Sleipner Skagerak
Rhoma 2.66 2.66 2.65 2.65 2.68
m 2 2 * * 2.02
Rw = 0.07 ohmm at 20 oC, Temp Gradient : 2.6  oC
Formation NetGros Phif Sw. Klogh_A * Klogh_H Klogh-G
Heather 0.109 0.118 0.521 0.001 0.001 0.001
Hugin 0.813 0.230 0.154 162 2.0 17
"""




def test_permeability_equation_coefficients():
    assert pe.parse_klogh_equation("10**(8* Phif-9*Vsh+2)") == (8.0, -9.0, 2.0)
    assert pe.parse_klogh_equation("10**(17.4*Phif-3*Vsh-1.85)") == (17.4, -3.0, -1.85)
    assert pe.parse_klogh_equation("something else") is None


def test_klogh_check_assigns_samples_to_the_reports_transform_and_counts_floor_and_misfits():
    eq = pd.DataFrame({"quantity": "permeability", "formation": ["Hugin", "Skagerak"], "equation": ["10**(8* Phif-9*Vsh+2)", "10**(17.4*Phif-3*Vsh-1.85)"]})
    phif = np.array([0.2, 0.25, 0.1, 0.15, 0.2])
    vsh = np.array([0.1, 0.2, 0.5, 0.3, 0.1])
    k = 10 ** (8 * phif - 9 * vsh + 2)
    k[2] = 0.001                                                                         # a sample at the floor
    k[3] = k[3] * 3                                                                      # a sample that follows no transform
    k[4] = 10 ** (17.4 * phif[4] - 3 * vsh[4] - 1.85)                                    # a Skagerrak sample
    out = pe.klogh_check(pd.DataFrame({"phif": phif, "vsh": vsh, "klogh": k}), eq, {"Skagerak": "Skagerrak"}, 0.05)
    n = out.set_index("transform")["samples"]
    assert n["Hugin"] == 2 and n["Skagerrak"] == 1 and n["floor value"] == 1 and n["none of the report's transforms"] == 1


def test_formation_intervals_use_aliases_merge_variants_and_end_at_td():
    tops = pd.DataFrame({"wellbore_id": 1, "strat_unit": ["Draupne Fm. Top", "Hugin Fm. Top", "Hugin Fm. Sand VOLVE Top", "Skagerak Fm. Top", "Other Top"], "top_md_m": [100.0, 120.0, 150.0, 180.0, 190.0]})
    iv = pe.formation_intervals(tops, 1, {"Skagerak": "Skagerrak"}, 250.0)
    assert iv["formation"].tolist() == ["Draupne", "Hugin", "Skagerrak"]
    assert iv["top_md_m"].tolist() == [100.0, 120.0, 180.0] and iv["base_md_m"].tolist() == [120.0, 180.0, 190.0]     # the Hugin variants merge; the next pick ends an interval


def test_zone_summary_counts_net_by_the_sand_flag_and_blanks_small_groups():
    md = np.arange(0.0, 100.0, 1.0)
    c = pd.DataFrame({"md_m": md, "phif": 0.2, "sw": 0.3, "klogh": 100.0, "sand_flag": (md % 2 == 0).astype(float)})
    z = pe.zone_summary(c, pd.DataFrame({"formation": ["Hugin", "Small"], "top_md_m": [0.0, 90.0], "base_md_m": [80.0, 92.0]}), min_n=5)
    h = z[z["formation"] == "Hugin"].iloc[0]
    assert abs(h["net_to_gross_flag"] - 0.5) < 1e-9 and abs(h["klogh_geometric_md"] - 100.0) < 1e-6 and not h["extends_to_last_sample"]
    assert np.isnan(z[z["formation"] == "Small"]["phif_net_mean"].iloc[0])                      # fewer than five net samples: blank, never zero
