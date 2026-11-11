"""Conventional core analysis: cleaning, depth shift, checks and summaries.  Pure functions; they return `(result, events)` where they change or flag data, and only the pipeline and its stage modules (`src/pipeline.py`, `src/stages/`) does file I/O.

What the data are: routine plug measurements of a laboratory (gas permeability `Kg`, the reciprocal mean pressure of the gas test `1/Pm`, the Klinkenberg-corrected permeability `Kl`, helium porosity, fluid saturations
by retort or extraction, grain density, a lithological description), for the horizontal and, in some samples, the vertical plug.  Depths are MD from the rotary table AS CORED; the core shift table gives, per core, the
top and base depth shifted to the log depth (company `UNKNOWN`, method not stated).  Rules: the laboratory's markers (`NMP`, `NPP`) are kept as flags and never read as numbers; no value is removed (outside a plausible range
it is flagged); both depths are kept and the shifted depth is labelled as interpolated from the shift table.  ASSUMPTIONS (listed, in the adapter): `1/Pm` is in 1/bar; the meaning of `NMP` and `NPP`.
"""
from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd

from ep_physics_base.events import dq_event

NUMERIC = ["kg_h_md", "inv_pm_h", "kl_h_md", "kg_v_md", "inv_pm_v", "kl_v_md", "por_h_pct", "por_v_pct", "por_sum_pct", "so_pct", "sw_pct", "grain_density_h_gcc", "grain_density_v_gcc"]
SCHEMA_VERSIONS = {"core_plug": "1.0", "core_shift": "1.0", "core_summary": "1.0", "core_klinkenberg": "1.0"}


def clean_conventional(raw: pd.DataFrame, cfg: dict[str, Any], well_name: str, source_file: str, sid: str) -> tuple[pd.DataFrame, list[dict]]:
    """One row per plug sample with typed numbers.  A laboratory marker (`cfg['markers']`: text -> meaning) in a measurement column becomes null and sets `<column>_marker`; a text cell that is neither a marker nor
    a number becomes null and is counted.  Blank cells are null.  The row order is the order of the file; nothing is dropped."""
    events: list[dict] = []
    d = raw.copy()
    markers = {str(k).strip().upper(): v for k, v in cfg.get("markers", {}).items()}
    for c in ("core_no", "sample_no", "depth_md_m"):
        d[c] = pd.to_numeric(d.get(c), errors="coerce")
    flagged: dict[str, int] = {}
    unreadable = 0
    for c in NUMERIC:
        if c not in d:
            d[c] = np.nan
            continue
        s = d[c].map(lambda v: v.strip() if isinstance(v, str) else v)
        mk = s.map(lambda v: markers.get(v.upper()) if isinstance(v, str) else None)
        d[f"{c}_marker"] = mk.astype("string")
        num = pd.to_numeric(s.where(mk.isna()), errors="coerce")
        unreadable += int((s.notna() & (s != "") & num.isna() & mk.isna()).sum())
        flagged[c] = int(mk.notna().sum())
        d[c] = num.astype("Float64")
    d["description"] = d.get("description", pd.Series(index=d.index, dtype=object)).map(lambda v: v.strip() if isinstance(v, str) and v.strip() else None).astype("string")
    d["well_name"] = well_name
    d["source_file"] = source_file
    d = d[d["depth_md_m"].notna()].reset_index(drop=True)
    d["core_no"] = d["core_no"].astype("Int64")
    d["sample_no"] = d["sample_no"].astype("Int64")
    n_mark = sum(flagged.values())
    if n_mark:
        events.append(dq_event(sid, "core_plug", "lab_marker_not_a_number", n_mark, "info", [f"{well_name}: {k} x{v}" for k, v in flagged.items() if v][:5],
                               "NMP / NPP written by the laboratory in a measurement cell: kept as a marker column, the value is null (meaning is an ASSUMPTION in the adapter)"))
    if unreadable:
        events.append(dq_event(sid, "core_plug", "text_not_a_number", unreadable, "warning", [well_name], "a text cell that is neither a marker nor a number was set to null"))
    return d, events


def qc_flags(d: pd.DataFrame, limits: dict[str, Any], sid: str) -> tuple[pd.DataFrame, list[dict]]:
    """Plausibility flags (limits are ASSUMPTIONS in the adapter); a flagged value stays in the table.  Columns: `flag_porosity_range`, `flag_grain_density`, `flag_saturation_sum`, `flag_klinkenberg_above_gas`
    (corrected permeability above the gas permeability, which a Klinkenberg correction cannot give), `flag_depth_order` (depth not increasing within a core)."""
    out = d.copy()
    p_lo, p_hi = limits["porosity_pct"]
    g_lo, g_hi = limits["grain_density_gcc"]
    por = out[["por_h_pct", "por_v_pct"]].astype("float64")
    out["flag_porosity_range"] = ((por < p_lo) | (por > p_hi)).any(axis=1)
    gd = out[["grain_density_h_gcc", "grain_density_v_gcc"]].astype("float64")
    out["flag_grain_density"] = ((gd < g_lo) | (gd > g_hi)).any(axis=1)
    out["flag_saturation_sum"] = (out["so_pct"].astype("float64") + out["sw_pct"].astype("float64")) > limits["saturation_sum_max_pct"]
    tol = limits.get("klinkenberg_tolerance", 1.02)
    out["flag_klinkenberg_above_gas"] = ((out["kl_h_md"].astype("float64") > tol * out["kg_h_md"].astype("float64")) | (out["kl_v_md"].astype("float64") > tol * out["kg_v_md"].astype("float64")))
    out["flag_depth_order"] = out.groupby(["well_name", "core_no"])["depth_md_m"].diff().le(0).fillna(False)
    events = []
    for col, what in [("flag_porosity_range", "porosity outside the plausible range"), ("flag_grain_density", "grain density outside the plausible range (heavy minerals such as siderite or pyrite, or an error)"),
                      ("flag_saturation_sum", "oil + water saturation above the limit"), ("flag_klinkenberg_above_gas", "Klinkenberg-corrected permeability above the gas permeability"), ("flag_depth_order", "depth not increasing within a core")]:
        n = int(out[col].sum())
        if n:
            ex = [f"{r.well_name} {r.depth_md_m:.2f}" for r in out[out[col]].head(5).itertuples()]
            events.append(dq_event(sid, "core_plug", col, n, "warning", ex, what + "; kept and flagged"))
    return out, events


def apply_core_shift(plug: pd.DataFrame, shift: pd.DataFrame) -> pd.DataFrame:
    """Adds `depth_shifted_md_m`, `shift_m` and `shift_basis`.  Per core the shift is interpolated linearly between the top and the base of the shift table (as cored -> shifted); a sample outside the cored interval of its core
    takes the nearest end shift (`shift_basis = 'end value'`); a sample whose core is not in the table, or a well without a table, keeps `depth_shifted_md_m` null (`shift_basis = 'no shift table'`).  The raw depth is never changed."""
    out = plug.copy()
    out["depth_shifted_md_m"] = np.nan
    out["shift_m"] = np.nan
    out["shift_basis"] = "no shift table"
    if shift is None or shift.empty:
        return out
    s = shift.copy()
    for c in ("core_no", "top_md_m", "base_md_m", "top_shifted_md_m", "base_shifted_md_m"):
        s[c] = pd.to_numeric(s[c], errors="coerce")
    for r in s.dropna(subset=["core_no", "top_md_m", "base_md_m", "top_shifted_md_m", "base_shifted_md_m"]).itertuples():
        m = (out["well_name"] == r.well_name) & (out["core_no"] == int(r.core_no))
        if not m.any():
            continue
        z = out.loc[m, "depth_md_m"].astype(float)
        top_sh, base_sh = r.top_shifted_md_m - r.top_md_m, r.base_shifted_md_m - r.base_md_m
        span = r.base_md_m - r.top_md_m
        w = ((z - r.top_md_m) / span).clip(0, 1) if span > 0 else pd.Series(0.0, index=z.index)
        sh = top_sh + w * (base_sh - top_sh)
        out.loc[m, "shift_m"] = sh
        out.loc[m, "depth_shifted_md_m"] = z + sh
        out.loc[m, "shift_basis"] = np.where((z < r.top_md_m) | (z > r.base_md_m), "end value", "interpolated")
    return out


def klinkenberg(plug: pd.DataFrame, pm_unit_per_bar: float = 1.0) -> pd.DataFrame:
    """The slip factor implied by each sample with gas permeability, corrected permeability and `1/Pm`: Kl = Kg / (1 + b / Pm), so b = (Kg / Kl - 1) x Pm [bar] with Pm = 1 / (1/Pm).  `pm_unit_per_bar` converts the
    unit of `1/Pm` if it is not 1/bar (ASSUMPTION: 1/bar).  Only samples with all three values and Kg above Kl; an implied b is a derived consistency check of the laboratory's correction, not a measurement."""
    rows = []
    for ori, kg, kl, ip in (("horizontal", "kg_h_md", "kl_h_md", "inv_pm_h"), ("vertical", "kg_v_md", "kl_v_md", "inv_pm_v")):
        g = plug[["well_name", "depth_md_m", kg, kl, ip]].dropna().astype({kg: float, kl: float, ip: float})
        g = g[(g[kg] > g[kl]) & (g[kl] > 0) & (g[ip] > 0)]
        pm = 1.0 / (g[ip] * pm_unit_per_bar)
        rows.append(pd.DataFrame({"well_name": g["well_name"], "depth_md_m": g["depth_md_m"], "orientation": ori, "kg_md": g[kg], "kl_md": g[kl], "inv_pm_per_bar": g[ip] * pm_unit_per_bar, "slip_b_bar": (g[kg] / g[kl] - 1.0) * pm}))
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def core_summary(plug: pd.DataFrame, min_n: int = 5) -> pd.DataFrame:
    """One row per well and core: samples, cored interval of the samples, porosity (mean, median), permeability (geometric mean and median of the horizontal corrected permeability) and the share of samples with a marker.  A
    statistic from fewer than `min_n` samples is null, never zero.  The geometric mean is the usual average of a permeability distribution that spans decades."""
    rows = []
    for (w, c), g in plug.groupby(["well_name", "core_no"], dropna=False):
        por = g["por_h_pct"].astype("float64").dropna()
        k = g["kl_h_md"].astype("float64").dropna()
        k = k[k > 0]
        rows.append({"well_name": w, "core_no": c, "samples": len(g), "top_md_m": g["depth_md_m"].astype(float).min(), "base_md_m": g["depth_md_m"].astype(float).max(),
                     "porosity_n": len(por), "porosity_mean_pct": por.mean() if len(por) >= min_n else np.nan, "porosity_median_pct": por.median() if len(por) >= min_n else np.nan,
                     "perm_n": len(k), "perm_geomean_md": float(np.exp(np.log(k).mean())) if len(k) >= min_n else np.nan, "perm_median_md": k.median() if len(k) >= min_n else np.nan,
                     "marker_share": float(g[[c_ for c_ in g.columns if c_.endswith("_marker")]].notna().any(axis=1).mean()), "flagged_share": float(g[[c_ for c_ in g.columns if c_.startswith("flag_")]].any(axis=1).mean())})
    return pd.DataFrame(rows)


def core_vs_log(plug: pd.DataFrame, curve: pd.DataFrame, max_gap_m: float, min_n: int = 5, perm_floor_md: float | None = None) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Plug values against the petrophysical interpretation (CPI) of the same well at the SHIFTED core depth.  Per plug with a horizontal porosity or permeability and a shifted depth, the log value is the nearest sample of the
    curve within `max_gap_m` (else null).  Returns the paired table (`phi_core`, `phi_log`, `k_core_md`, `k_log_md`, `sw_core`, `sw_log` as fractions or mD) and a summary per well (pairs, bias log minus core, RMSE, Pearson r for
    porosity; median ratio and Pearson r of log10 permeability).  The two are NOT independent: an interpretation is usually calibrated to the core of the same well, so agreement is a consistency check and a
    disagreement is a finding; a plug (cm) and a log (decimetres to metres) also sample different volumes.  A summary from fewer than `min_n` pairs is null.  `perm_floor_md`: an interpreted permeability that repeats exactly one value (a floor or cutoff of the interpretation, an inference from the data) is flagged (`k_log_at_floor`) and the ratio is also given for the pairs above it."""
    pairs = []
    sums = []
    for w, g in plug.groupby("well_name"):
        c = curve[curve["well_name"] == w].sort_values("md_m")
        g = g[g["depth_shifted_md_m"].notna()].copy()
        if c.empty or g.empty:
            continue
        idx = np.searchsorted(c["md_m"].to_numpy(), g["depth_shifted_md_m"].astype(float).to_numpy())
        idx = np.clip(idx, 1, len(c) - 1)
        left, right = c.iloc[idx - 1], c.iloc[idx]
        pick = np.where(np.abs(g["depth_shifted_md_m"].to_numpy(float) - left["md_m"].to_numpy()) <= np.abs(right["md_m"].to_numpy() - g["depth_shifted_md_m"].to_numpy(float)), idx - 1, idx)
        near = c.iloc[pick].reset_index(drop=True)
        gap = np.abs(near["md_m"].to_numpy() - g["depth_shifted_md_m"].to_numpy(float))
        ok = gap <= max_gap_m
        out = pd.DataFrame({"well_name": w, "depth_md_m": g["depth_md_m"].to_numpy(), "depth_shifted_md_m": g["depth_shifted_md_m"].to_numpy(float), "log_gap_m": gap,
                            "phi_core": g["por_h_pct"].astype("float64").to_numpy() / 100.0, "k_core_md": g["kl_h_md"].astype("float64").to_numpy(), "sw_core": g["sw_pct"].astype("float64").to_numpy() / 100.0})
        for col, name in (("phif", "phi_log"), ("klogh", "k_log_md"), ("sw", "sw_log"), ("vsh", "vsh_log")):
            out[name] = np.where(ok, near[col].to_numpy(float), np.nan) if col in near else np.nan
        out["k_log_at_floor"] = (out["k_log_md"] <= perm_floor_md * 1.0001) if perm_floor_md else False
        pairs.append(out)
        p = out.dropna(subset=["phi_core", "phi_log"])
        k = out.dropna(subset=["k_core_md", "k_log_md"])
        k = k[(k["k_core_md"] > 0) & (k["k_log_md"] > 0)]
        sums.append({"well_name": w, "phi_pairs": len(p), "phi_bias_log_minus_core": (p["phi_log"] - p["phi_core"]).mean() if len(p) >= min_n else np.nan,
                     "phi_rmse": float(np.sqrt(((p["phi_log"] - p["phi_core"]) ** 2).mean())) if len(p) >= min_n else np.nan, "phi_r": p["phi_log"].corr(p["phi_core"]) if len(p) >= min_n else np.nan,
                     "k_pairs": len(k), "k_median_ratio_log_over_core": float(np.median(k["k_log_md"] / k["k_core_md"])) if len(k) >= min_n else np.nan,
                     "log10k_r": np.log10(k["k_log_md"]).corr(np.log10(k["k_core_md"])) if len(k) >= min_n else np.nan,
                     "k_log_floor_share": float(k["k_log_at_floor"].mean()) if len(k) and perm_floor_md else np.nan})
        ka = k[~k["k_log_at_floor"]] if perm_floor_md else k
        sums[-1]["k_pairs_above_floor"] = len(ka)
        sums[-1]["k_median_ratio_above_floor"] = float(np.median(ka["k_log_md"] / ka["k_core_md"])) if len(ka) >= min_n else np.nan
        sums[-1]["log10k_r_above_floor"] = np.log10(ka["k_log_md"]).corr(np.log10(ka["k_core_md"])) if len(ka) >= min_n else np.nan
    return (pd.concat(pairs, ignore_index=True) if pairs else pd.DataFrame()), pd.DataFrame(sums)


def lab_klinkenberg(kn2: pd.Series, inv_pm: pd.Series, f: dict[str, float]) -> pd.Series:
    """The laboratory's own empirical Klinkenberg correction (final report of 15/9-19 A, PDF page 11, image checked): below `low_max` mD `Kl = low_a x Kn2 ** low_b`; above it two iterations of
    `m = m_a x Kl ** m_b`, `Kl = Kn2 - m / Pm` starting from the low-range value, with `1/Pm` as delivered.  Returns the permeability the formula gives for each plug (null where an input is null); compared with the
    delivered `Kl` this reproduces the laboratory's arithmetic and tests the unit of `1/Pm`."""
    kn2 = kn2.astype("float64")
    ip = inv_pm.astype("float64")
    kl = f["low_a"] * kn2 ** f["low_b"]
    for _ in range(int(f.get("iterations", 2))):
        kl_hi = kn2 - f["m_a"] * kl ** f["m_b"] * ip
        kl = kl_hi.where(kn2 > f["low_max"], f["low_a"] * kn2 ** f["low_b"])
    return kl.where(kn2.notna() & ip.notna())


def klinkenberg_reproduction(plug: pd.DataFrame, f: dict[str, float], tol: float = 0.02) -> pd.DataFrame:
    """Per plug orientation: the delivered corrected permeability, the laboratory formula's value and their ratio; `agrees` is true within `tol`.  Only plugs with all three values."""
    rows = []
    for ori, kg, ip, kl in (("horizontal", "kg_h_md", "inv_pm_h", "kl_h_md"), ("vertical", "kg_v_md", "inv_pm_v", "kl_v_md")):
        g = plug[["well_name", "depth_md_m", kg, ip, kl]].dropna().copy()
        calc = lab_klinkenberg(g[kg], g[ip], f)
        rows.append(pd.DataFrame({"well_name": g["well_name"], "depth_md_m": g["depth_md_m"], "orientation": ori, "kn2_md": g[kg].astype(float), "inv_pm": g[ip].astype(float), "kl_delivered_md": g[kl].astype(float),
                                  "kl_formula_md": calc, "ratio_formula_over_delivered": calc / g[kl].astype(float), "agrees": (calc / g[kl].astype(float) - 1).abs() <= tol}))
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
