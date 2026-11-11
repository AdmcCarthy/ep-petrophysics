
import numpy as np
import pandas as pd

from ep_petrophysics import core_analysis as ca

CFG = {"markers": {"NMP": "no measurement possible (assumed)", "NPP": "no plug possible (assumed)"}}
LIM = {"porosity_pct": [0.0, 40.0], "grain_density_gcc": [2.55, 2.95], "saturation_sum_max_pct": 100.0, "klinkenberg_tolerance": 1.02}


def _raw():
    return pd.DataFrame({
        "core_no": [1, 1, 1, 2], "sample_no": [1, 2, 3, 4], "depth_md_m": [100.0, 100.5, 101.0, 110.0],
        "kg_h_md": ["13.8", "NMP ", "  ", "50"], "inv_pm_h": [0.66, None, None, 0.9], "kl_h_md": ["11.5", "NMP ", "  ", "60"],
        "por_h_pct": [17.0, 14.8, "NPP ", 50.0], "grain_density_h_gcc": [2.66, 3.4, 2.65, 2.65], "so_pct": [None, None, None, 60.0], "sw_pct": [None, None, None, 55.0],
        "description": ["Sst", " ", None, "Sst"],
    })


def test_markers_are_flags_not_numbers_and_nothing_is_dropped():
    d, ev = ca.clean_conventional(_raw(), CFG, "W", "f.xls", "s")
    assert len(d) == 4
    assert d["kg_h_md"].isna().tolist() == [False, True, True, False]            # NMP and blanks are null
    assert d["kg_h_md_marker"].notna().tolist() == [False, True, False, False]   # only the marker is recorded
    assert d["por_h_pct_marker"].iloc[2] == "no plug possible (assumed)"
    assert any(e["rule"] == "lab_marker_not_a_number" for e in ev)


def test_qc_flags_keep_the_values():
    d, _ = ca.clean_conventional(_raw(), CFG, "W", "f.xls", "s")
    q, ev = ca.qc_flags(d, LIM, "s")
    assert q["flag_grain_density"].tolist() == [False, True, False, False]
    assert q["flag_porosity_range"].tolist() == [False, False, False, True]         # 50 % flagged, still present
    assert q["flag_saturation_sum"].tolist() == [False, False, False, True]
    assert q["flag_klinkenberg_above_gas"].tolist() == [False, False, False, True]  # Kl 60 above Kg 50
    assert q["por_h_pct"].iloc[3] == 50.0


def test_core_shift_interpolates_and_labels_end_values():
    plug = pd.DataFrame({"well_name": ["W"] * 4, "core_no": [1, 1, 1, 2], "depth_md_m": [100.0, 105.0, 111.0, 120.0]})
    shift = pd.DataFrame({"well_name": ["W"], "core_no": [1], "top_md_m": [100.0], "base_md_m": [110.0], "top_shifted_md_m": [101.0], "base_shifted_md_m": [113.0]})
    out = ca.apply_core_shift(plug, shift)
    assert out["shift_m"].round(2).tolist()[:3] == [1.0, 2.0, 3.0]                  # 1 m at the top, 3 m at the base, 111 m is past the base: end value
    assert out["shift_basis"].tolist() == ["interpolated", "interpolated", "end value", "no shift table"]
    assert np.isnan(out["depth_shifted_md_m"].iloc[3])                                # a core without a row is not shifted
    assert out["depth_md_m"].tolist() == [100.0, 105.0, 111.0, 120.0]                 # the raw depth is never changed


def test_klinkenberg_slip_factor_round_trip():
    # Kl = Kg / (1 + b / Pm): b = 0.5 bar, Pm = 2 bar -> Kl = Kg / 1.25
    plug = pd.DataFrame({"well_name": ["W"], "depth_md_m": [1.0], "kg_h_md": [10.0], "kl_h_md": [8.0], "inv_pm_h": [0.5], "kg_v_md": [np.nan], "kl_v_md": [np.nan], "inv_pm_v": [np.nan]})
    k = ca.klinkenberg(plug)
    assert len(k) == 1 and abs(k["slip_b_bar"].iloc[0] - 0.5) < 1e-9


def test_core_vs_log_pairs_at_shifted_depth_and_flags_the_floor():
    plug = pd.DataFrame({"well_name": ["W"] * 3, "depth_md_m": [1.0, 2.0, 3.0], "depth_shifted_md_m": [1.1, 2.1, 3.1], "por_h_pct": [20.0, 15.0, 10.0], "kl_h_md": [100.0, 10.0, 1.0], "sw_pct": [np.nan] * 3})
    curve = pd.DataFrame({"well_name": ["W"] * 3, "md_m": [1.1, 2.1, 3.1], "phif": [0.19, 0.16, 0.1], "klogh": [90.0, 0.1, 0.1], "sw": [0.2] * 3, "vsh": [0.1] * 3})
    pairs, summ = ca.core_vs_log(plug, curve, 0.2, min_n=2, perm_floor_md=0.1)
    assert pairs["phi_log"].round(2).tolist() == [0.19, 0.16, 0.1]
    assert pairs["k_log_at_floor"].tolist() == [False, True, True]
    assert summ["k_pairs_above_floor"].iloc[0] == 1





