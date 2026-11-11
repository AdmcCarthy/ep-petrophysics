import numpy as np
import pandas as pd

from ep_petrophysics import facies as fa

BY_KEY = {"159F14": 5351}


def _raw():
    return pd.DataFrame({"uwi": "NO 15/9-F-14", "well_name_file": "15/9-F-14", "curve_type": ["GENETIC FACIES"] * 4, "source": "GM09", "top_md_m": [100.0, 110.0, 125.0, 130.0], "base_md_m": [110.0, 120.0, 130.0, 140.0],
                         "facies_class": ["F-TIDAL BAR", "F-TIDAL CHANNEL", "F-TIDAL BAR", "F-MARSH"], "rock_percent": 100.0})


def test_clean_facies_links_the_well_and_sees_the_gap():
    d, ev = fa.clean_facies(_raw(), BY_KEY, "s")
    assert (d["wellbore_id"] == 5351).all() and d["is_genetic_facies"].all()
    assert d["gap_above_m"].isna().tolist() == [True, True, False, True] and d.loc[2, "gap_above_m"] == 5.0      # 120 -> 125 leaves a 5 m gap


def test_gap_and_overlap_are_flagged_correctly():
    raw = _raw()
    raw.loc[2, "top_md_m"] = 118.0                       # overlaps the previous base (120) by 2 m
    raw.loc[3, "top_md_m"] = 133.0                       # leaves a 3 m gap after the previous base (130)
    raw.loc[2, "base_md_m"] = 130.0
    d, ev = fa.clean_facies(raw, BY_KEY, "s")
    assert d.loc[2, "overlap_above_m"] == 2.0 and d.loc[3, "gap_above_m"] == 3.0
    assert {e["rule"] for e in ev} == {"gap_above_m", "overlap_above_m"}


def test_assign_interval_uses_top_inclusive_base_exclusive_and_leaves_outside_empty():
    lab = fa.assign_interval(np.array([99.0, 100.0, 109.9, 110.0, 140.0]), np.array([100.0, 110.0]), np.array([110.0, 120.0]), np.array(["A", "B"], dtype=object))
    assert lab.tolist() == [None, "A", "A", "B", None]


def test_core_by_facies_assigns_only_inside_the_interval_and_blanks_small_groups():
    g, _ = fa.clean_facies(_raw(), BY_KEY, "s")
    g["well_key"] = "159F14"
    plug = pd.DataFrame({"well_name": "15/9-F-14", "depth_md_m": np.arange(100.0, 160.0, 1.0)})
    plug["depth_shifted_md_m"] = plug["depth_md_m"]
    plug["por_h_pct"] = 20.0
    plug["kl_h_md"] = np.where(plug["depth_md_m"].between(110, 119.99), 1000.0, 10.0)
    p, s = fa.core_by_facies(plug, g, min_n=5)
    assert p["facies"].notna().sum() == 35 and p.loc[p["depth_md_m"] == 150.0, "facies"].iloc[0] is None          # 100-120 and 125-140 are inside intervals, 120-125 is a gap
    ch = s[s["facies"] == "F-TIDAL CHANNEL"].iloc[0]
    assert ch["plugs"] == 10 and abs(ch["perm_geomean_md"] - 1000.0) < 1e-6
    big = fa.core_by_facies(plug.iloc[:104], g, min_n=20)[1]
    assert big["perm_geomean_md"].isna().all()                                                                      # fewer than min_n plugs: blank, never zero


def test_dip_depth_unit_is_decided_by_the_survey_not_the_label():
    md = np.arange(0.0, 4000.0, 50.0)
    survey = pd.DataFrame({"md_m": md, "inclination_deg": np.clip(md / 50.0, 0, 60)})
    dips = pd.DataFrame({"md_stated": [1000.0, 2000.0, 3000.0]})
    dips["deviation_deg"] = np.interp(dips["md_stated"], survey["md_m"], survey["inclination_deg"])      # the numbers are metres
    r = fa.dip_depth_unit_check(dips, survey)
    assert r["unit_of_numbers"] == "metres" and r["metres"] < 1e-9 and r["feet"] > 1


def test_circular_mean_wraps_around_north():
    assert abs(fa.circular_mean_deg(pd.Series([350.0, 10.0])) - 0.0) < 1e-6 or abs(fa.circular_mean_deg(pd.Series([350.0, 10.0])) - 360.0) < 1e-6
    assert abs(fa.circular_mean_deg(pd.Series([90.0, 100.0])) - 95.0) < 1e-6
