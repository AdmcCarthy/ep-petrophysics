"""Cased-hole logs of the Volve wells: cement bond (integrity logs) and production logging (pressure, temperature, fluid density, spinner).  Pure functions; only the pipeline and its stage modules (`src/pipeline.py`, `src/stages/`) reads and writes.

Cement: the provider's own indices (CBL bond index `BI`, ultrasonic bond index `USBI`, gas index `USGI`, average acoustic impedance `AIAV`) are summarised per run; the classes (good, moderate, poor) are ASSUMED limits in
`config/borehole.yaml` (`cement`), a bond log is NOT a barrier verdict (no standard acceptance is applied: requests req.101 and req.102).  Production logging: the measured pressure gradient of the wellbore fluid on the date of a
run replaces an assumed gradient where one exists, and the log pressure is compared with the permanent gauge on the same day.
"""
from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd

SCHEMA_VERSIONS = {"well_cement_profile": "1.0", "well_cement_summary": "1.0", "well_plt_profile": "1.0", "well_plt_gradient": "1.0", "well_plt_vs_gauge": "1.0"}


def _one_channel(curves: pd.DataFrame, role: str, priority: list[str]) -> pd.DataFrame:
    """For one member and role the channel with most bins (ties: the order of the adapter's mnemonic list); columns md_m, value."""
    g = curves[curves["role"] == role]
    if g.empty:
        return pd.DataFrame(columns=["md_m", "value"])
    n = g.groupby("channel").size()
    order = {m.upper(): i for i, m in enumerate(priority)}
    best = sorted(n.index, key=lambda c: (-n[c], order.get(str(c).upper(), 99)))[0]
    return g[g["channel"] == best].groupby("md_m", as_index=False)["value"].median()


def cement_profile(curves: pd.DataFrame, spec: dict[str, Any]) -> pd.DataFrame:
    """One wide table per member on a common 1 m grid: a column per cement role (the channel with most bins), `log_date`, `well_folder`.  A run that reads several frames or repeat passes is merged by the median per bin."""
    out = []
    for member, g in curves.groupby("member"):
        wide = None
        for role, r in spec.items():
            if not isinstance(r, dict) or "mnemonics" not in r:
                continue
            c = _one_channel(g, role, r["mnemonics"])
            if c.empty:
                continue
            c = c.rename(columns={"value": role})
            wide = c if wide is None else wide.merge(c, on="md_m", how="outer")
        if wide is not None:
            wide["member"], wide["log_date"], wide["well_folder"] = member, g["log_date"].iloc[0], g["well_folder"].iloc[0]
            out.append(wide.sort_values("md_m"))
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame()


def _longest_run(md: np.ndarray, ok: np.ndarray, gap_m: float) -> tuple[float, float, float]:
    """The longest stretch (m) of `ok` bins in which neighbouring ok bins are at most `gap_m` apart; (length, top, base); the bins are 1 m, a stretch of one bin counts as its bin width."""
    idx = np.flatnonzero(ok)
    if idx.size == 0:
        return 0.0, np.nan, np.nan
    best = (0.0, np.nan, np.nan)
    start = prev = idx[0]
    for i in list(idx[1:]) + [None]:
        if i is None or md[i] - md[prev] > gap_m:
            length = md[prev] - md[start] + 1.0
            if length > best[0]:
                best = (float(length), float(md[start]), float(md[prev]))
            if i is not None:
                start = i
        if i is not None:
            prev = i
    return best


def cement_summary(profile: pd.DataFrame, cfg: dict[str, Any]) -> pd.DataFrame:
    """Per run: logged interval, which indices exist, the shares of the logged bins in each bond class of the ultrasonic index (else the CBL index; classes and the allowed gap are ASSUMED limits), the longest stretch of good bond
    and where it lies, the median impedance and eccentering.  `index_used` names the index; a run without a bond index is listed with `index_used = none`."""
    good, moderate, gap = float(cfg["good_bond"]), float(cfg["moderate_bond"]), float(cfg["gap_m"])
    rows = []
    for member, g in profile.groupby("member"):
        idx_col = next((c for c in ("bond_index_usit", "bond_index_cbl") if c in g and g[c].notna().sum() > 0), None)
        r = {"member": member, "well_folder": g["well_folder"].iloc[0], "log_date": g["log_date"].iloc[0], "top_md_m": g["md_m"].min(), "base_md_m": g["md_m"].max(), "bins": len(g), "index_used": idx_col or "none"}
        for c in ("bond_index_cbl", "bond_index_usit", "gas_index_usit", "impedance", "cbl_amplitude", "eccentering"):
            r[f"{c}_bins"] = int(g[c].notna().sum()) if c in g else 0
        if idx_col:
            v = g[idx_col]
            ok = v.notna()
            n = int(ok.sum())
            top_share = float(v[ok].round(3).value_counts(normalize=True).iloc[0]) if n else np.nan                                  # the share of bins that carry the single most common value
            r["index_constant_share"] = top_share
            r["index_looks_constant"] = bool(top_share > float(cfg.get("constant_share", 0.8)))                                      # a bond index that does not vary is a default or a channel that was not computed, not a measurement
            r.update({"good_share": float((v[ok] >= good).mean()), "moderate_share": float(((v[ok] >= moderate) & (v[ok] < good)).mean()), "poor_share": float((v[ok] < moderate).mean()), "index_median": float(v[ok].median()), "indexed_bins": n})
            length, top, base = _longest_run(g["md_m"].to_numpy(float), (v >= good).fillna(False).to_numpy(), gap)
            r.update({"longest_good_m": length, "longest_good_top_md_m": top, "longest_good_base_md_m": base})
        for c, k in (("impedance", "impedance_median_mrayl"), ("eccentering", "eccentering_median_in"), ("gas_index_usit", "gas_index_median")):
            r[k] = float(g[c].median()) if c in g and g[c].notna().any() else np.nan
        rows.append(r)
    out = pd.DataFrame(rows)
    out["limitations"] = ("provider indices as delivered (their processing and flags are not interpreted); the classes are assumed limits, not an acceptance criterion; a bond log shows the casing-cement interface, not zonal isolation; an index that is the same value in most bins (`index_looks_constant`) is not a measurement; "
                          "cement type, slurry and the date of the cement job are not in the data read")
    return out


def plt_profile(curves: pd.DataFrame, spec: dict[str, Any], tvdss_at) -> pd.DataFrame:
    """Per well folder and log date the pressure, temperature and fluid density on a 5 m grid pooled over the passes of the day (the median and the spread between passes), with the TVD below sea level from the selected survey
    (`tvdss_at(well_folder, md_array)`, NaN outside the surveyed interval).  The pressure used is the first available of the roles `pressure_well`, `pressure_borehole`, `pressure_quartz` (named in `pressure_role`)."""
    out = []
    for (folder, date), g in curves.groupby(["well_folder", "log_date"]):
        for prole in ("pressure_well", "pressure_borehole", "pressure_quartz"):
            p = g[g["role"] == prole]
            if len(p):
                break
        else:
            continue
        p = p.assign(md5=(np.floor(p["md_m"] / 5.0) * 5.0 + 2.5))
        grid = p.groupby("md5").agg(pressure_bar=("value", "median"), pressure_spread_bar=("value", lambda s: float(s.max() - s.min())), passes=("member", "nunique")).reset_index().rename(columns={"md5": "md_m"})
        for role, name in (("temperature", "temperature_c"), ("fluid_density", "fluid_density_g_cm3"), ("spinner", "spinner_c_s")):
            q = g[g["role"] == role]
            if len(q):
                q = q.assign(md5=(np.floor(q["md_m"] / 5.0) * 5.0 + 2.5)).groupby("md5")["value"].median().rename(name).reset_index().rename(columns={"md5": "md_m"})
                grid = grid.merge(q, on="md_m", how="left")
        grid["well_folder"], grid["log_date"], grid["pressure_role"] = folder, date, prole
        grid["tvdss_m"] = np.asarray(tvdss_at(folder, grid["md_m"].to_numpy(float)), dtype=float)
        out.append(grid)
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame()


def _fluid_class(slope: float, cfg: dict[str, Any]) -> str:
    """The class of a measured pressure gradient by the equivalent density it implies, with the configured limits (ASSUMPTIONS: gas below `gas_below`, water-like from `water_from` to `water_to`, above that not a static fluid gradient)."""
    if not np.isfinite(slope):
        return "not fitted"
    if slope < float(cfg["gas_below_bar_per_m"]):
        return "gas-like (below the oil range)"
    if slope < float(cfg["water_from_bar_per_m"]):
        return "oil or gas-oil mixture"
    if slope <= float(cfg["water_to_bar_per_m"]):
        return "water-like"
    return "above a water gradient: not a static fluid column (pressure changed in time, or a heavy brine)"


def plt_gradient(prof: pd.DataFrame, cfg: dict[str, Any]) -> pd.DataFrame:
    """The measured pressure gradient of the wellbore fluid in windows of `window_tvd_m` of TVD (step half a window) per well folder and log date: a robust straight line (Theil-Sen) of pressure against TVD below sea level over the bins
    of the window that have a TVD and a pressure above `min_pressure_bar`.  One line over a whole well would average gas, oil and water, so the windows are the result.  A window is `fitted` only when it has `min_bins` bins, at
    least `min_tvd_span_m` of vertical span and a mean pressure spread between passes below `max_spread_bar` (when the passes disagree the pressure changed in time and the profile is not static).  The equivalent density is gradient /
    0.0980665 [sg]; `fluid_class` classes it; `fluid_density_median_g_cm3` is the log's own density channel for the same window (a different measurement, kept for comparison)."""
    from scipy.stats import theilslopes

    rows = []
    w = float(cfg["window_tvd_m"])
    for (folder, date), g in prof.groupby(["well_folder", "log_date"]):
        d = g[g["tvdss_m"].notna() & g["pressure_bar"].gt(float(cfg["min_pressure_bar"]))].sort_values("tvdss_m")
        if d.empty:
            continue
        starts = np.arange(d["tvdss_m"].min(), max(d["tvdss_m"].max() - w / 2, d["tvdss_m"].min() + 1.0), w / 2)
        for z0 in starts:
            s = d[(d["tvdss_m"] >= z0) & (d["tvdss_m"] < z0 + w)]
            span = float(s["tvdss_m"].max() - s["tvdss_m"].min()) if len(s) else 0.0
            r = {"well_folder": folder, "log_date": date, "tvdss_top_m": float(s["tvdss_m"].min()) if len(s) else np.nan, "tvdss_base_m": float(s["tvdss_m"].max()) if len(s) else np.nan, "md_top_m": float(s["md_m"].min()) if len(s) else np.nan,
                 "md_base_m": float(s["md_m"].max()) if len(s) else np.nan, "bins": len(s), "tvd_span_m": span, "pressure_role": g["pressure_role"].iloc[0], "passes_max": int(s["passes"].max()) if len(s) else 0,
                 "mean_pass_spread_bar": float(s["pressure_spread_bar"].mean()) if len(s) else np.nan,
                 "fluid_density_median_g_cm3": float(s["fluid_density_g_cm3"].median()) if len(s) and "fluid_density_g_cm3" in s and s["fluid_density_g_cm3"].notna().any() else np.nan}
            if len(s) < int(cfg["min_bins"]) or span < float(cfg["min_tvd_span_m"]):
                r.update({"status": "no gradient: vertical span or bins too small"})
            elif r["mean_pass_spread_bar"] > float(cfg["max_spread_bar"]):
                r.update({"status": "no gradient: the passes disagree (pressure changed in time)"})
            else:
                slope, icpt, lo, hi = theilslopes(s["pressure_bar"].to_numpy(float), s["tvdss_m"].to_numpy(float))
                r.update({"gradient_bar_per_m": float(slope), "gradient_low_bar_per_m": float(lo), "gradient_high_bar_per_m": float(hi), "equivalent_density_sg": float(slope / 0.0980665), "status": "fitted"})
            r["fluid_class"] = _fluid_class(r.get("gradient_bar_per_m", np.nan), cfg)
            rows.append(r)
    out = pd.DataFrame(rows)
    for c in ("gradient_bar_per_m", "gradient_low_bar_per_m", "gradient_high_bar_per_m", "equivalent_density_sg"):
        if c not in out:
            out[c] = np.nan
    return out


def plt_vs_gauge(prof: pd.DataFrame, gauge_daily: pd.DataFrame, gauges: pd.DataFrame, folder_wellbore: dict[str, int], max_distance_m: float, gradients: pd.DataFrame | None = None, spread_max_bar: float = 3.0, min_passes: int = 2) -> pd.DataFrame:
    """The production-log pressure near the permanent gauge against the gauge's daily average on the same day.  Only a STATIONARY pressure counts: a bin read by at least `min_passes` passes of the day whose spread is at most
    `spread_max_bar` (a bin read once, while the tool moves, can carry surge or a different time of the build-up; such a pair is listed with that reason, never compared).  The bin nearest to the gauge MD within `max_distance_m` is
    carried to the gauge TVD with the measured gradient of the nearest fitted window of the same day (`log_pressure_at_gauge_tvd_bar`, derived; blank without a gradient).  The hours on stream of that day tell whether the gauge
    read a flowing or a shut-in pressure.  Both numbers are a day's values, not one instant."""
    rows = []
    gd = gauges.set_index("wellbore_id") if len(gauges) else pd.DataFrame()
    for (folder, date), g in prof.groupby(["well_folder", "log_date"]):
        wid = folder_wellbore.get(folder)
        if wid is None or wid not in gd.index:
            continue
        gmd, gz = float(gd.loc[wid, "gauge_md_m"]), gd.loc[wid, "gauge_tvdss_m"]
        day = gauge_daily[(gauge_daily["wellbore_id"] == wid) & (gauge_daily["date"] == pd.Timestamp(date))]
        gauge_bar = float(day["downhole_pressure_bar"].iloc[0]) if len(day) and day["downhole_pressure_bar"].notna().any() else np.nan
        r = {"well_folder": folder, "log_date": date, "wellbore_id": wid, "gauge_md_m": gmd, "gauge_daily_bar": gauge_bar, "hours_on_stream": float(day["on_stream_hrs"].iloc[0]) if len(day) else np.nan}
        stat = g[(g["passes"] >= min_passes) & (g["pressure_spread_bar"] <= spread_max_bar)]
        near_any = (g["md_m"] - gmd).abs().min()
        if stat.empty or (stat["md_m"] - gmd).abs().min() > max_distance_m:
            r.update({"status": "no stationary log pressure near the gauge" + (" (the log reaches the gauge depth on a single moving pass only)" if near_any <= max_distance_m else " (the log does not reach the gauge depth)")})
            rows.append(r)
            continue
        i = (stat["md_m"] - gmd).abs().idxmin()
        z_log = stat.loc[i, "tvdss_m"]
        r.update({"log_bin_md_m": float(stat.loc[i, "md_m"]), "distance_m": float(abs(stat.loc[i, "md_m"] - gmd)), "log_pressure_bar": float(stat.loc[i, "pressure_bar"]), "log_passes": int(stat.loc[i, "passes"]), "log_spread_bar": float(stat.loc[i, "pressure_spread_bar"])})
        grad = gradients[(gradients["well_folder"] == folder) & (gradients["log_date"] == date) & gradients["gradient_bar_per_m"].notna()] if gradients is not None and len(gradients) else pd.DataFrame()
        if len(grad):
            grad = grad.assign(_d=((grad["tvdss_top_m"] + grad["tvdss_base_m"]) / 2 - gz).abs()).sort_values("_d").head(1)
        r["gradient_used_bar_per_m"] = float(grad["gradient_bar_per_m"].iloc[0]) if len(grad) else np.nan
        r["log_pressure_at_gauge_tvd_bar"] = r["log_pressure_bar"] - r["gradient_used_bar_per_m"] * (z_log - gz) if pd.notna(r["gradient_used_bar_per_m"]) and pd.notna(z_log) and pd.notna(gz) else np.nan
        r["difference_at_gauge_tvd_bar"] = r["log_pressure_at_gauge_tvd_bar"] - gauge_bar
        r["status"] = "compared" if pd.notna(gauge_bar) else "no gauge value that day"
        rows.append(r)
    return pd.DataFrame(rows)