"""Petrophysical interpretation outputs: formation intervals, zone summaries, reproduction of the report's permeability transform and comparison with the report's zone table.  Pure functions.

What the data are: the interpretation of the operator's petrophysicists (porosity from density with a neutron correction calibrated to overburden-corrected core, clay volume from gamma ray, Archie water saturation, a permeability
transform regressed on core): output logs (`PHIF`, `SW`, `VSH`, `KLOGH`, `BVW` and the flags `SAND_FLAG`, `CARB_FLAG`, `COAL_FLAG`) and, for 9 wellbores, the model parameters and the zone results of the report.  Rules: the interpretation is a
MODEL calibrated to core, never a measurement; the report's equations are applied to the output logs as a CHECK of the arithmetic (a disagreement is a finding); the net criterion of the report is not stated, so the zone summary uses the
interpretation's sand flag as net (an ASSUMPTION) and is compared with the report's table, never substituted for it; a formation interval is the interval between interpretation picks (the deepest extends to the last curve sample).
"""
from __future__ import annotations

import re
from typing import Any

import numpy as np
import pandas as pd

from ep_physics_base.names import name_key, unit_key_of_name

SCHEMA_VERSIONS = {"petro_curve": "1.0", "petro_model": "1.0", "petro_zone": "1.0", "petro_check": "1.0"}
FORMATIONS = ("DRAUPNE", "HEATHER", "HUGIN", "SLEIPNER", "SKAGERRAK")


def formation_intervals(tops: pd.DataFrame, wellbore_id: int, aliases: dict[str, str], td_md_m: float) -> pd.DataFrame:
    """The formation intervals of a wellbore from its picks (`strat_unit`, `top_md_m`): each pick named like one of the report formations (after the aliases, `Skagerak` = `Skagerrak`) runs from its top to the next pick's top
    (the deepest to `td_md_m`).  Picks that are not one of the formations still end the previous interval.  Returns `formation`, `top_md_m`, `base_md_m`."""
    t = tops[tops["wellbore_id"] == int(wellbore_id)].sort_values("top_md_m")
    t = t.drop_duplicates(["strat_unit", "top_md_m"])
    if t.empty:
        return pd.DataFrame(columns=["formation", "top_md_m", "base_md_m"])
    alias = {name_key(k): name_key(v) for k, v in aliases.items()}
    rows = []
    tops_m = t["top_md_m"].to_numpy(float)
    for i, r in enumerate(t.itertuples()):
        key = unit_key_of_name(r.strat_unit)
        key = alias.get(key, key)
        if key in FORMATIONS:
            rows.append({"formation": key.title(), "top_md_m": float(r.top_md_m), "base_md_m": float(tops_m[i + 1]) if i + 1 < len(tops_m) else float(td_md_m)})
    d = pd.DataFrame(rows)
    if d.empty:
        return d
    merged = [d.iloc[0].to_dict()]
    for r in d.iloc[1:].to_dict("records"):                                   # consecutive picks of one formation (variants such as `Hugin Fm. Sand`) form one interval
        if r["formation"] == merged[-1]["formation"] and abs(r["top_md_m"] - merged[-1]["base_md_m"]) < 0.01:
            merged[-1]["base_md_m"] = r["base_md_m"]
        else:
            merged.append(r)
    return pd.DataFrame(merged)


def _harmonic(v: pd.Series) -> float:
    v = v.dropna()
    v = v[v > 0]
    return float(len(v) / (1.0 / v).sum()) if len(v) else np.nan


def zone_summary(curve: pd.DataFrame, intervals: pd.DataFrame, min_n: int = 5) -> pd.DataFrame:
    """Per formation interval: samples, the share of samples flagged sand (net to gross as the flag gives it) and, over the flagged samples, the mean porosity and water saturation and the arithmetic, harmonic and geometric mean
    permeability.  A sample needs a porosity to count; a statistic from fewer than `min_n` flagged samples is null.  Where the curve has no sand flag the flagged statistics are null."""
    rows = []
    for r in intervals.itertuples():
        g = curve[(curve["md_m"] >= r.top_md_m) & (curve["md_m"] < r.base_md_m) & curve["phif"].notna()]
        row = {"formation": r.formation, "top_md_m": r.top_md_m, "base_md_m": r.base_md_m, "samples": len(g), "extends_to_last_sample": bool(r.base_md_m >= float(curve["md_m"].max()) - 0.2)}
        if "sand_flag" in g and g["sand_flag"].notna().any() and len(g):
            s = g[g["sand_flag"] == 1]
            row.update({"net_to_gross_flag": float(len(s) / len(g)), "net_samples": len(s)})
            ok = len(s) >= min_n
            row["phif_net_mean"] = float(s["phif"].mean()) if ok else np.nan
            row["sw_net_mean"] = float(s["sw"].mean()) if ok else np.nan
            k = s["klogh"].dropna()
            k = k[k > 0]
            row["klogh_arithmetic_md"] = float(k.mean()) if ok and len(k) else np.nan
            row["klogh_harmonic_md"] = _harmonic(k) if ok else np.nan
            row["klogh_geometric_md"] = float(np.exp(np.log(k).mean())) if ok and len(k) else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def compare_with_report(zone: pd.DataFrame, results: pd.DataFrame) -> pd.DataFrame:
    """Each report zone row beside the zone summary computed from the logs (same formation name): differences and relative differences of porosity, water saturation, net to gross and the three permeability means."""
    if zone.empty or results.empty:
        return pd.DataFrame()
    z = zone.copy()
    z["formation_key"] = z["formation"].map(name_key)
    r = results.copy()
    r["formation_key"] = r["formation"].map(lambda f: {"SKAGERAK": "SKAGERRAK"}.get(name_key(f), name_key(f)))
    m = r.merge(z, on="formation_key", how="inner", suffixes=("_report", "_logs"))
    out = pd.DataFrame({"formation": m["formation_report"], "net_to_gross_report": m["net_to_gross"], "net_to_gross_logs": m["net_to_gross_flag"], "phif_report": m["phif"], "phif_logs": m["phif_net_mean"],
                        "sw_report": m["sw"], "sw_logs": m["sw_net_mean"], "klogh_arithmetic_report_md": m["klogh_arithmetic_md_report"], "klogh_arithmetic_logs_md": m["klogh_arithmetic_md_logs"],
                        "klogh_harmonic_report_md": m["klogh_harmonic_md_report"], "klogh_harmonic_logs_md": m["klogh_harmonic_md_logs"], "klogh_geometric_report_md": m["klogh_geometric_md_report"],
                        "klogh_geometric_logs_md": m["klogh_geometric_md_logs"], "samples": m["samples"], "extends_to_last_sample": m["extends_to_last_sample"]})
    out["phif_diff"] = out["phif_logs"] - out["phif_report"]
    out["sw_diff"] = out["sw_logs"] - out["sw_report"]
    out["net_to_gross_diff"] = out["net_to_gross_logs"] - out["net_to_gross_report"]
    out["klogh_geometric_ratio_logs_over_report"] = out["klogh_geometric_logs_md"] / out["klogh_geometric_report_md"]
    return out


_EQ = re.compile(r"10\*\*\(\s*([-+]?\d+(?:\.\d+)?)\s*\*\s*Phif\s*([-+]\s*\d+(?:\.\d+)?)\s*\*\s*Vsh\s*([-+]\s*\d+(?:\.\d+)?)\s*\)", re.I)


def parse_klogh_equation(text: str) -> tuple[float, float, float] | None:
    """The coefficients (a, b, c) of a transform `10**(a*Phif + b*Vsh + c)` written in the report (`10**(8* Phif-9*Vsh+2)`); None when the text has another form."""
    m = _EQ.search(text)
    if not m:
        return None
    return float(m.group(1)), float(m.group(2).replace(" ", "")), float(m.group(3).replace(" ", ""))


def klogh_check(curve: pd.DataFrame, equations: pd.DataFrame, aliases: dict[str, str], tol: float = 0.05, floor_md: float = 0.001) -> pd.DataFrame:
    """Which of the report's permeability transforms produced each output permeability sample.  Every transform of the report (one per formation: `10**(a*Phif + b*Vsh + c)`) is applied to the output porosity and clay volume;
    a sample belongs to the first transform that reproduces its output permeability within `tol` (relative); a sample at the floor value `floor_md` is counted as floor; the rest as `none`.  Formation intervals are not used (the picks
    do not cover every formation boundary): the check asks whether the output permeability is the report's arithmetic, sample by sample.  Returns one row per transform plus `floor` and `none` with the number and share of samples."""
    eq = []
    for e in equations[equations["quantity"] == "permeability"].itertuples():
        co = parse_klogh_equation(e.equation)
        f = alias_key(e.formation, aliases)
        if co and f and f not in [x[0] for x in eq]:
            eq.append((f, co))
    g = curve.dropna(subset=["phif", "vsh", "klogh"])
    g = g[g["phif"].notna()]
    n = len(g)
    rows = []
    if not n or not eq:
        return pd.DataFrame()
    taken = np.zeros(n, dtype=bool)
    for f, co in eq:
        ratio = 10 ** (co[0] * g["phif"].to_numpy() + co[1] * g["vsh"].to_numpy() + co[2]) / np.where(g["klogh"].to_numpy() > 0, g["klogh"].to_numpy(), np.nan)
        hit = (np.abs(ratio - 1) <= tol) & ~taken
        taken |= hit
        rows.append({"transform": f.title(), "coefficients": f"{co[0]:g} x PHIF {co[1]:+g} x VSH {co[2]:+g}", "samples": int(hit.sum()), "share": float(hit.sum() / n)})
    fl = (np.abs(g["klogh"].to_numpy() - floor_md) <= 1e-9) & ~taken
    taken |= fl
    rows.append({"transform": "floor value", "coefficients": f"{floor_md:g} mD", "samples": int(fl.sum()), "share": float(fl.sum() / n)})
    rows.append({"transform": "none of the report's transforms", "coefficients": "", "samples": int((~taken).sum()), "share": float((~taken).sum() / n)})
    return pd.DataFrame(rows)


def alias_key(name: Any, aliases: dict[str, str]) -> str:
    k = name_key(name)
    a = {name_key(x): name_key(y) for x, y in aliases.items()}
    return a.get(k, k)
