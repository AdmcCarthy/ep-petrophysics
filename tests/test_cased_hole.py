"""Cement bond summaries and production-logging pressure gradients (synthetic data only)."""
import numpy as np
import pandas as pd

from ep_petrophysics import cased_hole as ch

CEM = {"bond_index_usit": {"mnemonics": ["USBI", "USBI_RF"]}, "bond_index_cbl": {"mnemonics": ["BI"]}, "impedance": {"mnemonics": ["AIAV"]}}
CFG = {"good_bond": 0.8, "moderate_bond": 0.5, "gap_m": 2.0}


def _curves(values, channel="USBI", role="bond_index_usit"):
    md = np.arange(len(values), dtype=float) + 1000.5
    return pd.DataFrame({"member": "m1", "well_folder": "W", "log_date": "2013-06-05", "channel": channel, "role": role, "md_m": md, "value": values, "n": 5})


def test_cement_summary_classes_and_longest_good_stretch_with_a_gap():
    v = [0.9] * 10 + [0.3] + [0.9] * 5 + [0.2] * 4 + [0.6] * 5                    # good 10 m, one poor bin, good 5 m, poor 4 m, moderate 5 m
    c = pd.concat([_curves(v), _curves([3.0] * len(v), "AIAV", "impedance")], ignore_index=True)
    s = ch.cement_summary(ch.cement_profile(c, CEM), CFG).iloc[0]
    assert s["index_used"] == "bond_index_usit" and s["bins"] == 25
    assert abs(s["good_share"] - 15 / 25) < 1e-9 and abs(s["moderate_share"] - 5 / 25) < 1e-9 and abs(s["poor_share"] - 5 / 25) < 1e-9
    assert s["longest_good_m"] == 16.0 and s["impedance_median_mrayl"] == 3.0            # the 1 m poor bin (gap 2 m) does not break the stretch: 10 + 1 + 5
    assert "not a barrier" in s["limitations"] or "not zonal isolation" in s["limitations"] or "zonal isolation" in s["limitations"]


def test_cement_without_a_bond_index_is_listed_not_dropped():
    c = _curves([3.0] * 5, "AIAV", "impedance")
    s = ch.cement_summary(ch.cement_profile(c, CEM), CFG).iloc[0]
    assert s["index_used"] == "none" and s["impedance_median_mrayl"] == 3.0


PCFG = {"min_pressure_bar": 50.0, "min_bins": 8, "min_tvd_span_m": 50.0, "window_tvd_m": 200.0, "max_spread_bar": 3.0, "gas_below_bar_per_m": 0.04, "water_from_bar_per_m": 0.085, "water_to_bar_per_m": 0.115}


def _plt(p0, grad, n=40, step=5.0, vertical=True):
    md = 3000.0 + step * np.arange(n)
    z = (md - 55.0) if vertical else 2900.0 + 0.1 * np.arange(n)
    p = p0 + grad * (z - z[0])
    rows = [{"member": "m1", "well_folder": "W", "log_date": "2014-11-28", "channel": "WPRE", "role": "pressure_well", "md_m": m + 0.5, "value": v, "n": 3} for m, v in zip(md, p)]
    rows += [{"member": "m1", "well_folder": "W", "log_date": "2014-11-28", "channel": "MWFD", "role": "fluid_density", "md_m": m + 0.5, "value": 0.55, "n": 3} for m in md]
    rows += [{**r, "member": "m2", "value": r["value"] + 0.1} for r in rows]                          # a second pass: the bins are read twice (stationary)
    return pd.DataFrame(rows), z


def test_plt_gradient_recovers_a_known_gradient_and_refuses_a_horizontal_section():
    c, z = _plt(250.0, 0.062)
    prof = ch.plt_profile(c, {}, lambda f, md: np.asarray(md) - 55.5)
    gw = ch.plt_gradient(prof, PCFG)
    g = gw[gw["status"] == "fitted"].iloc[0]
    assert g["fluid_class"].startswith("oil") and g["status"] == "fitted" and abs(g["gradient_bar_per_m"] - 0.062) < 1e-3 and abs(g["equivalent_density_sg"] - 0.062 / 0.0980665) < 0.02 and abs(g["fluid_density_median_g_cm3"] - 0.6) < 1e-9
    prof2 = prof.assign(tvdss_m=2900.0 + 0.1 * np.arange(len(prof)))                        # almost no vertical span
    assert ch.plt_gradient(prof2, PCFG).iloc[0]["status"].startswith("no gradient")


def test_plt_vs_gauge_difference_and_distance_limit():
    c, z = _plt(250.0, 0.062)
    prof = ch.plt_profile(c, {}, lambda f, md: np.asarray(md) - 55.5)
    grad = ch.plt_gradient(prof, PCFG)
    gauges = pd.DataFrame({"wellbore_id": [1], "gauge_md_m": [2950.0], "gauge_tvdss_m": [2950.0 - 55.5]})
    daily = pd.DataFrame({"wellbore_id": pd.array([1], dtype="Int64"), "date": [pd.Timestamp("2014-11-28")], "downhole_pressure_bar": [247.0], "on_stream_hrs": [0.0]})
    r = ch.plt_vs_gauge(prof, daily, gauges, {"W": 1}, 150.0, grad).iloc[0]
    assert r["status"] == "compared" and abs(r["distance_m"] - 52.5) < 1e-6 and r["hours_on_stream"] == 0.0 and r["log_passes"] == 2
    assert abs(r["log_pressure_at_gauge_tvd_bar"] - (250.0 + 0.062 * (2894.5 - 2945.0))) < 0.2            # carried to the gauge depth with the measured gradient
    far = ch.plt_vs_gauge(prof, daily, gauges.assign(gauge_md_m=2000.0), {"W": 1}, 150.0, grad).iloc[0]
    assert far["status"].startswith("no stationary log pressure") and "does not reach" in far["status"]
    single = c[c["member"] == "m1"]                                                                         # one pass only: never compared
    prof1 = ch.plt_profile(single, {}, lambda f, md: np.asarray(md) - 55.5)
    assert "single moving pass" in ch.plt_vs_gauge(prof1, daily, gauges, {"W": 1}, 150.0, None).iloc[0]["status"]
