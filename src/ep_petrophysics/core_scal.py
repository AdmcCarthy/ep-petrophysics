"""Special core analysis (SCAL): end points, relative permeability curves, wettability and the checks against the routine core data.  Pure functions.

What the data are: a laboratory's end-point table (initial water saturation `Swi`, residual oil saturation `Sor`, oil permeability at `Swi`, water permeability at `Sor`, end-point water relative permeability, Amott
wettability indices), the oil relative permeability curves of centrifuge imbibition tests (ambient pressure, 70 degC, crude oil), and the unsteady-state water floods of composite cores at reservoir conditions.  Rules: the
laboratory's values are shown as delivered; recomputed values (`krw(Sor)`, the Corey exponent) are labelled derived and are checks, never replacements; a disagreement is a finding; no value is removed.
ASSUMPTIONS (adapter): the depth `mRKB` of the overview is the depth of the core as cored (MD from the rotary table), the same system as the routine plug table; a three-composite test is labelled by the depth of
its plugs rounded to the metre.
"""
from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd

from ep_physics_base.events import dq_event

SCHEMA_VERSIONS = {"core_scal_endpoint": "1.0", "core_relperm_curve": "1.0", "core_scal_check": "1.0"}


def family_of_test(test_type: Any) -> str:
    """The family of a test from the laboratory's text (rules: `USS` or `Res.cond` = water flood at reservoir conditions; `Cent` = centrifuge relative permeability; `USBM` = wettability test with its plug state)."""
    t = str(test_type).lower()
    if "uss" in t or "res.cond" in t:
        return "water flood, reservoir conditions"
    if "cent" in t:
        return "centrifuge, imbibition oil relative permeability"
    if "usbm" in t:
        return "wettability, " + ("fresh" if "fresh" in t else "restored" if "restored" in t else "cleaned" if "cleaned" in t else "state not stated")
    return "not classified"


def clean_endpoints(raw: pd.DataFrame, well_name: str, source_file: str, sid: str, tol: float = 0.03) -> tuple[pd.DataFrame, list[dict]]:
    """Adds the family, the recomputed end-point water relative permeability `krw_sor_check = kw(Sor) / ko(Swi)` and the flag `flag_krw_differs` when it differs from the stated value by more than `tol`.  The Amott
    indices must lie between 0 and 1 (`flag_index_range`).  Nothing is removed."""
    d = raw.copy()
    d["well_name"] = well_name
    d["source_file"] = source_file
    d["test_family"] = d["test_type"].map(family_of_test)
    d["krw_sor_check"] = d["kw_sor_md"] / d["ko_swi_md"]
    d["flag_krw_differs"] = (d["krw_sor_check"] - d["krw_sor"]).abs() > tol
    d["flag_index_range"] = ((d["wwi"] < 0) | (d["wwi"] > 1) | (d["owi"] < 0) | (d["owi"] > 1)).fillna(False)
    d["wi_amott"] = d["wwi"] - d["owi"]                                  # the laboratory's own definition (report p53): WIAmott = WWI - OWI
    events = []
    n = int(d["flag_krw_differs"].sum())
    if n:
        events.append(dq_event(sid, "core_scal", "krw_sor_differs_from_kw_over_ko", n, "warning", [f"{r.depth_md_m:.2f}" for r in d[d["flag_krw_differs"]].head(5).itertuples()],
                               "end-point krw(Sor) differs from kw(Sor) / ko(Swi) by more than the tolerance; both are kept"))
    return d, events


def waterflood_checks(curves: pd.DataFrame, endpoints: pd.DataFrame, sid: str, tol: float = 0.02) -> tuple[pd.DataFrame, list[dict]]:
    """Per water-flood curve (by the nominal depth): the first water saturation of the curve against the stated `Swi` of the end-point table (they must agree), the last saturation against `1 - Sor`, and the number of
    points.  A curve whose starting saturation is not the stated `Swi` is flagged (`flag_swi_mismatch`): the sheet may carry a copied starting point."""
    rows, events = [], []
    ep = endpoints[endpoints["test_family"] == "water flood, reservoir conditions"]
    for depth, g in curves.groupby("nominal_depth_md_m"):
        e = ep[(ep["depth_md_m"] - depth).abs() < 0.5]
        if e.empty:
            continue
        e = e.iloc[0]
        s = g.dropna(subset=["sw"])
        first, last = float(g["sw"].dropna().iloc[0]), float(s["sw"].iloc[-1])
        rows.append({"nominal_depth_md_m": depth, "points": len(g), "sw_first": first, "swi_stated": e["swi"], "sw_last": last, "one_minus_sor_stated": 1 - e["sor"],
                     "flag_swi_mismatch": abs(first - e["swi"]) > tol, "flag_sor_mismatch": abs(last - (1 - e["sor"])) > 0.05})
    out = pd.DataFrame(rows)
    for col, what in (("flag_swi_mismatch", "the first saturation of the curve is not the stated initial water saturation"), ("flag_sor_mismatch", "the last saturation of the curve is not 1 - the stated residual oil saturation")):
        if len(out) and out[col].any():
            events.append(dq_event(sid, "core_scal", col, int(out[col].sum()), "warning", [f"{r:.0f}" for r in out.loc[out[col], "nominal_depth_md_m"]], what + "; kept"))
    return out, events


def fit_corey_oil(curves: pd.DataFrame, endpoints: pd.DataFrame, swi_from: str = "swi", sor_from: str = "sor") -> pd.DataFrame:
    """Corey exponent of the oil relative permeability of each centrifuge test, fitted by least squares on `kro = ((1 - Sw - Sor) / (1 - Swi - Sor)) ** n` over the points with `kro` between 0.001 and 0.99, with `Swi` and
    `Sor` of the end-point table (the report's own Sor of the curve is its last `Sw`: `Sor = 1 - Sw_max`).  Compared with the laboratory's stated exponent.  Derived; a fit is not a measurement."""
    rows = []
    for depth, g in curves.groupby("depth_md_m"):
        e = endpoints[(endpoints["depth_md_m"] - depth).abs() < 0.005]
        if e.empty:
            continue
        e = e.iloc[0]
        swi = float(g["sw"].min())
        sor = float(1 - g["sw"].max())
        x = (1 - g["sw"] - sor) / (1 - swi - sor)
        m = (g["kro"] > 0.001) & (g["kro"] < 0.99) & (x > 0.001) & (x < 0.999)
        if m.sum() < 8:
            continue
        n = float(np.sum(np.log(x[m]) * np.log(g.loc[m, "kro"])) / np.sum(np.log(x[m]) ** 2))        # through-the-origin fit of log kro = n log x
        pred = x[m] ** n
        rows.append({"depth_md_m": depth, "points": int(m.sum()), "corey_oil_fit": n, "corey_oil_stated": e["corey_oil"], "swi_curve": swi, "sor_curve": sor,
                     "swi_stated": e["swi"], "sor_stated": e["sor"], "rmse_kro": float(np.sqrt(np.mean((pred - g.loc[m, "kro"]) ** 2)))})
    return pd.DataFrame(rows)


def scal_vs_routine(endpoints: pd.DataFrame, plug: pd.DataFrame, window_m: float = 0.15) -> pd.DataFrame:
    """Each SCAL plug against the routine plug nearest to it in as-cored depth (within `window_m`): the laboratory's `kL` and helium porosity of the SCAL plug against the routine Klinkenberg-corrected permeability
    and porosity.  The SCAL plugs were cleaned again and measured at the end of the SCAL program, the routine plugs at the first analysis: a difference is a finding about the depth match, the repeatability and the
    cleaning, not an error of either.  The routine depth is the as-cored depth (the same system as the overview, an assumption)."""
    rows = []
    for r in endpoints.dropna(subset=["kl_md"]).itertuples():
        g = plug[(plug["well_name"] == r.well_name) & ((plug["depth_md_m"].astype(float) - r.depth_md_m).abs() <= window_m)]
        if g.empty:
            continue
        j = (g["depth_md_m"].astype(float) - r.depth_md_m).abs().idxmin()
        p = g.loc[j]
        rows.append({"well_name": r.well_name, "depth_md_m": r.depth_md_m, "test_family": r.test_family, "routine_depth_md_m": float(p["depth_md_m"]), "gap_m": float(abs(p["depth_md_m"] - r.depth_md_m)),
                     "kl_scal_md": r.kl_md, "kl_routine_md": float(p["kl_h_md"]) if pd.notna(p["kl_h_md"]) else np.nan, "phi_scal": r.porosity_frac, "phi_routine": float(p["por_h_pct"]) / 100 if pd.notna(p["por_h_pct"]) else np.nan})
    d = pd.DataFrame(rows)
    if len(d):
        d["kl_ratio_scal_over_routine"] = d["kl_scal_md"] / d["kl_routine_md"]
        d["phi_diff_scal_minus_routine"] = d["phi_scal"] - d["phi_routine"]
    return d
