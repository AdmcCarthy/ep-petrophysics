"""Rock-physics logs (LFP): physical and arithmetic checks, the effect of the pore fluid on the elastic properties, and the interpreter's assumptions against other data.  Pure functions.

What the data are: the elastic properties (velocities, impedances, Poisson's ratio, density) of the logs and their fluid-substituted versions for brine (`_B`), oil (`_O`) and gas (`_G`) of the operator's rock-physics workflow, with the
petrophysical inputs and the interpreter's listing of parameters (2002).  Rules: the curves are shown as delivered; the substituted curves are a MODEL (Gassmann-type substitution, assumptions in the primer), not measurements; physical checks
(impedance = density x velocity, velocity from the sonic, the ordering brine >= oil >= gas of the velocity and impedance in a sand) are tests of the product, and a failure is a finding, never corrected; the sand flag stands for net sand (assumption);
implausible values are flagged and counted, not removed.
"""
from __future__ import annotations

import re
from typing import Any

import numpy as np
import pandas as pd

SCHEMA_VERSIONS = {"lfp_curve": "1.0", "lfp_document": "1.0", "lfp_check": "1.0", "lfp_fluid_effect": "1.0", "lfp_report_check": "1.0"}


def physical_checks(curve: pd.DataFrame, cfg: dict[str, Any]) -> pd.DataFrame:
    """Per wellbore name: samples with elastic properties and, over them, (1) the share where `LFP_AI` is within the tolerance of `LFP_RHOB x LFP_VP` (g/cm3 x m/s gives kPa.s/m), (2) the share where `LFP_VP` is within the tolerance
    of `304800 / LFP_DT` (sonic slowness in us/ft), (3) the share where `LFP_VPVS` equals `LFP_VP / LFP_VS`, (4) the share of samples flagged bad data, and (5) the counts of implausible values (impedance below the limit, effective
    porosity or gamma ray above theirs).  Columns that a wellbore lacks give null shares."""
    chk = cfg["checks"]
    rows = []
    for w, g in curve.groupby("well_name_file"):
        row = {"well_name_file": w, "samples": len(g)}
        ai = g["lfp_ai"].astype("float64")
        rho, vp = g["lfp_rhob"].astype("float64"), g["lfp_vp"].astype("float64")
        ok = ai.notna() & rho.notna() & vp.notna() & (rho > 0) & (vp > 0)
        row["impedance_samples"] = int(ok.sum())
        row["impedance_within_tol"] = float(((ai[ok] / (rho[ok] * vp[ok]) - 1).abs() <= chk["impedance_rel"]).mean()) if ok.any() else np.nan
        dt = g["lfp_dt"].astype("float64")
        ok2 = vp.notna() & dt.notna() & (dt > 0) & (vp > 0)
        row["vp_from_dt_samples"] = int(ok2.sum())
        row["vp_from_dt_within_tol"] = float(((vp[ok2] / (304800.0 / dt[ok2]) - 1).abs() <= chk["vp_from_dt_rel"]).mean()) if ok2.any() else np.nan
        vs, r = g["lfp_vs"].astype("float64"), g["lfp_vpvs"].astype("float64")
        ok3 = vp.notna() & vs.notna() & r.notna() & (vs > 0)
        row["vpvs_within_tol"] = float(((r[ok3] / (vp[ok3] / vs[ok3]) - 1).abs() <= 0.01).mean()) if ok3.any() else np.nan
        bad = g["lfp_baddata"].astype("float64")
        row["baddata_share"] = float((bad == 1).sum() / len(g)) if bad.notna().any() else np.nan                  # the flag is set (1) where the data are bad and left empty or 0 elsewhere
        row["ai_below_limit"] = int((ai < chk["plausible"]["ai_min"]).sum())
        row["phie_above_limit"] = int((g["lfp_phie"].astype("float64") > chk["plausible"]["phie_max"]).sum()) if "lfp_phie" in g else 0
        row["gr_above_limit"] = int((g["lfp_gr"].astype("float64") > chk["plausible"]["gr_max"]).sum()) if "lfp_gr" in g else 0
        rows.append(row)
    return pd.DataFrame(rows)


def fluid_ordering(curve: pd.DataFrame, intervals: dict[str, tuple[float, float]], min_n: int = 20) -> pd.DataFrame:
    """In each wellbore's net sand (sand flag 1, bad data 0, porosity above 0.05) of the stated interval: the share of samples where the substituted velocity and impedance follow the order brine >= oil and oil >= gas, separately
    (Gassmann substitution raises the fluid bulk modulus and density of brine above oil above gas).  `intervals` maps a wellbore name to the (top, base) of the Hugin interval from the interpreter's listing."""
    rows = []
    for w, g in curve.groupby("well_name_file"):
        top, base = intervals.get(w, (np.nan, np.nan))
        s = g[(g["md_m"] >= top) & (g["md_m"] < base)] if np.isfinite(top) else g.iloc[0:0]
        if "lfp_sand" in s:
            s = s[(s["lfp_sand"] == 1)]
        if "lfp_baddata" in s:
            s = s[~(s["lfp_baddata"] == 1)]
        s = s[s["lfp_phie"].astype("float64") > 0.05] if "lfp_phie" in s else s
        row = {"well_name_file": w, "interval_top_md_m": top, "interval_base_md_m": base, "net_sand_samples": len(s)}
        for q in ("vp", "ai"):
            b, o, gg = (s[f"lfp_{q}_{x}"].astype("float64") for x in ("b", "o", "g"))
            ok = b.notna() & o.notna() & gg.notna()
            row[f"{q}_brine_ge_oil"] = float((b[ok] >= o[ok] - 1e-6).mean()) if ok.sum() >= min_n else np.nan
            row[f"{q}_oil_ge_gas"] = float((o[ok] >= gg[ok] - 1e-6).mean()) if ok.sum() >= min_n else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def fluid_effect(curve: pd.DataFrame, intervals: dict[str, tuple[float, float]], min_n: int = 20) -> pd.DataFrame:
    """The effect of replacing brine by oil or gas in the net sand of the interval, per wellbore: the median relative change of P-wave velocity, density, P-impedance and Vp/Vs between the substituted cases
    (`oil / brine - 1`, `gas / brine - 1`) and the number of samples.  A statistic from fewer than `min_n` samples is null.  Model results, not measurements."""
    rows = []
    for w, g in curve.groupby("well_name_file"):
        top, base = intervals.get(w, (np.nan, np.nan))
        s = g[(g["md_m"] >= top) & (g["md_m"] < base)] if np.isfinite(top) else g.iloc[0:0]
        if "lfp_sand" in s:
            s = s[s["lfp_sand"] == 1]
        if "lfp_baddata" in s:
            s = s[~(s["lfp_baddata"] == 1)]
        s = s[s["lfp_phie"].astype("float64") > 0.05] if "lfp_phie" in s else s
        for prop, label in (("vp", "P-wave velocity"), ("rhob", "density"), ("ai", "P-impedance"), ("vpvs", "Vp/Vs")):
            b = s[f"lfp_{prop}_b"].astype("float64")
            for case, col in (("oil", "o"), ("gas", "g")):
                x = s[f"lfp_{prop}_{col}"].astype("float64")
                ok = b.notna() & x.notna() & (b > 0)
                rows.append({"well_name_file": w, "property": label, "case": case, "samples": int(ok.sum()), "median_change_vs_brine_pct": float(100 * (x[ok] / b[ok] - 1).median()) if ok.sum() >= min_n else np.nan,
                             "p10_pct": float(100 * (x[ok] / b[ok] - 1).quantile(0.1)) if ok.sum() >= min_n else np.nan, "p90_pct": float(100 * (x[ok] / b[ok] - 1).quantile(0.9)) if ok.sum() >= min_n else np.nan})
    return pd.DataFrame(rows)


def hugin_interval(units: pd.DataFrame, name: str = "Hugin Fm.") -> tuple[float, float]:
    """The depth range of a listed unit (the interpreter's listing of units and depths)."""
    u = units[units["unit"].str.strip().str.lower() == name.lower()]
    return (float(u["top_md_m"].iloc[0]), float(u["base_md_m"].iloc[0])) if len(u) else (np.nan, np.nan)


def report_checks(curve: pd.DataFrame, fact: pd.DataFrame, params: pd.DataFrame, cfg: dict[str, Any]) -> pd.DataFrame:
    """What the Word LFP report states against the curves of the LAS file, per wellbore that has both (`cfg['report_to_curve']`: report well hint -> curve well name).  Checks: the evaluation interval stated against the curve depth range; the oil interval stated (`LFP_OIL 1:a-b m MD`)
    against the depths where the model flag is 3 (hydrocarbon); the sand flag against the stated limits (porosity above the sand porosity limit and shale volume below the sand shale limit; the porosity curve is not named in the report, so effective and total porosity are both tried);
    `LFP_SUBST = LFP_SAND`; the temperature and pressure straight lines (intercept + slope x depth) against the curves, tested with the depth as MD RKB and as MD less the stated RKB-MSL (the depth basis is not stated, so both residuals are kept and the better is reported).
    One row per check: `stated`, `found`, `agreement` (share or median residual with unit), `status` (`agrees` within tolerance or `differs`).  A difference is a finding, never a correction; the curves are the interpreter's model."""
    tol = cfg["report_check"]
    rows = []

    def add(w, check, stated, found, value, unit, status, note=""):
        rows.append({"well_name": w, "check": check, "stated": stated, "found": found, "agreement": value, "unit": unit, "status": status, "note": note})

    for hint, cname in cfg["report_to_curve"].items():
        f = fact[fact["well_name_file"] == hint]
        c = curve[curve["well_name_file"] == cname].sort_values("md_m")
        if f.empty or c.empty:
            continue
        f = f.iloc[0]
        md = c["md_m"].to_numpy(float)
        if pd.notna(f["evaluation_start_md_m"]) or pd.notna(f["evaluation_end_md_m"]):
            s, e = f["evaluation_start_md_m"], f["evaluation_end_md_m"]
            d = max(abs(md.min() - s) if pd.notna(s) else 0.0, abs(md.max() - e) if pd.notna(e) else 0.0)
            add(cname, "evaluation_interval", f"{s if pd.notna(s) else 'not stated'} to {e if pd.notna(e) else 'not stated'} m MD RKB", f"{md.min():.1f} to {md.max():.1f} m", d, "m (largest end difference)", "agrees" if d <= tol["interval_m"] else "differs")
        pr = params[(params["well_name_file"] == hint) & params["parameter"].eq("LFP_OIL")]
        if len(pr) and "lfp_model" in c:
            m = re.search(r"1:\s*(\d+(?:\.\d+)?)\s*-\s*(\d+(?:\.\d+)?)", str(pr["value_text"].iloc[0]))
            hc = c.loc[c["lfp_model"] == 3, "md_m"]
            if m:
                a, b = float(m.group(1)), float(m.group(2))
                d = max(abs(hc.min() - a), abs(hc.max() - b)) if len(hc) else np.inf
                add(cname, "oil_interval_vs_model_flag", f"{a:g} to {b:g} m MD", f"{hc.min():.1f} to {hc.max():.1f} m" if len(hc) else "no hydrocarbon samples", d if np.isfinite(d) else np.nan, "m (largest end difference)", "agrees" if d <= tol["interval_m"] else "differs", "model flag 3 = hydrocarbon (LFP_MODEL = 3 - 2 x LFP_WATER)")
            else:
                add(cname, "oil_interval_vs_model_flag", str(pr["value_text"].iloc[0]), f"{int((c['lfp_model'] == 3).sum())} hydrocarbon samples", float((c["lfp_model"] == 3).sum()), "samples", "agrees" if (c["lfp_model"] == 3).sum() == 0 else "differs", "no oil interval stated; the flag is 3 only where hydrocarbon")
        if {"lfp_sand", "lfp_vsh"} <= set(c.columns):
            for pc in ("lfp_phie", "lfp_phit"):
                ok = c[["lfp_sand", pc, "lfp_vsh"]].dropna()
                pred = ((ok[pc] > f["sand_phi_limit"]) & (ok["lfp_vsh"] < f["sand_vsh_limit"])).astype(int)
                add(cname, f"sand_flag_rule_{pc[4:]}", f"sand = porosity > {f['sand_phi_limit']:g} and Vsh < {f['sand_vsh_limit']:g}", f"{len(ok)} samples", float((pred == ok["lfp_sand"].astype(int)).mean()), "share of samples that follow the rule", "agrees" if (pred == ok["lfp_sand"].astype(int)).mean() >= tol["rule_share"] else "differs", "the report does not name the porosity curve of the rule")
        if {"lfp_subst", "lfp_sand"} <= set(c.columns):
            ok = c[["lfp_subst", "lfp_sand"]].dropna()
            add(cname, "subst_equals_sand", "LFP_SUBST = LFP_SAND", f"{len(ok)} samples", float((ok["lfp_subst"] == ok["lfp_sand"]).mean()), "share of samples equal", "agrees" if (ok["lfp_subst"] == ok["lfp_sand"]).mean() >= tol["rule_share"] else "differs")
        for var, icol, scol, unit in (("temperature", "temperature_intercept", "temperature_slope_per_m", "degC"), ("pressure", "pressure_intercept", "pressure_slope_per_m", "bar")):
            col = "lfp_temp" if var == "temperature" else "lfp_press"
            if col not in c or pd.isna(f[icol]):
                continue
            ok = c[["md_m", col]].dropna()
            res = {name: np.median(np.abs(ok[col].to_numpy() - (f[icol] + f[scol] * (ok["md_m"].to_numpy() - off)))) for name, off in (("MD RKB", 0.0), ("MD less RKB-MSL", float(f["rkb_msl_m"]) if pd.notna(f["rkb_msl_m"]) else np.nan)) if not np.isnan(off)}
            best = min(res, key=res.get)
            add(cname, f"{var}_line", f"{f[icol]:.4g} + {f[scol]:.5g} x depth", f"{len(ok)} samples", float(res[best]), f"{unit} (median |curve - line|, depth = {best})", "agrees" if res[best] <= tol[f"{var}_{unit}"] else "differs", "; ".join(f"{k}: {v:.3g}" for k, v in res.items()))
    return pd.DataFrame(rows, columns=["well_name", "check", "stated", "found", "agreement", "unit", "status", "note"])