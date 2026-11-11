"""Interpretation products: genetic facies and lithology intervals of the geomodel, the permeability logs, core by facies and the borehole-image dip interpretation.  Pure functions; they return `(result, events)` where they change or flag data.

What the data are: facies intervals interpreted by the geomodelling team (source `GM09`, the 2009 geomodel: genetic facies of the Hugin Formation interval; source `SFINX`: a lithology class per fine interval of the whole well), a permeability log made for the
geomodel (`KLOGH_NEW`, from the petrophysical interpretation) and a dip interpretation of an azimuthal image (bed boundaries, possible faults, unconformable surfaces) of 15/9-F-14.  Rules: intervals are shown as delivered (a gap or an overlap between
two intervals is flagged, not closed); a facies is interpretation, not measurement; statistics by facies need at least `min_n` samples and are blank otherwise; a plug is assigned to a facies only through its SHIFTED depth and only inside the
interpreted interval; the depth unit of the dip workbook is tested against the survey, not trusted from its label; dip azimuths are averaged as angles (circular mean).
"""
from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd

from ep_physics_base.names import name_key
from ep_physics_base.events import dq_event

SCHEMA_VERSIONS = {"facies_interval": "1.0", "facies_stats": "1.0", "core_by_facies": "1.0", "image_dip": "1.0", "image_dip_bin": "1.0"}


def clean_facies(raw: pd.DataFrame, by_key: dict[str, int], sid: str) -> tuple[pd.DataFrame, list[dict]]:
    """Typed intervals with the wellbore link (exact normalised name or the Volve name map, `link_basis`), the thickness, and the flags `gap_above_m` (a gap to the previous interval of the same curve) and `overlap_above_m`."""
    d = raw.copy()
    d["well_name_file"] = d["well_name_file"].astype(str).str.strip()
    d["wellbore_id"] = d["well_name_file"].map(lambda n: by_key.get(name_key(n))).astype("Int64")
    d["link_basis"] = np.where(d["wellbore_id"].notna(), "name key (exact or Volve name map)", "no Sodir wellbore of this name")
    d["thickness_m"] = d["base_md_m"] - d["top_md_m"]
    d = d.sort_values(["well_name_file", "curve_type", "top_md_m"]).reset_index(drop=True)
    prev_base = d.groupby(["well_name_file", "curve_type"])["base_md_m"].shift()
    d["gap_above_m"] = (d["top_md_m"] - prev_base).where((d["top_md_m"] - prev_base) > 0.01)
    d["overlap_above_m"] = (prev_base - d["top_md_m"]).where((prev_base - d["top_md_m"]) > 0.01)
    d["is_genetic_facies"] = d["curve_type"].astype(str).str.upper().str.contains("GENETIC")
    events = []
    for col, what in (("gap_above_m", "gap between two intervals of one curve"), ("overlap_above_m", "overlap between two intervals of one curve")):
        n = int(d[col].notna().sum())
        if n:
            events.append(dq_event(sid, "facies_interval", col, n, "info", sorted(d.loc[d[col].notna(), "well_name_file"].unique())[:5], what + ": kept as delivered"))
    return d, events


def assign_interval(depth: np.ndarray, top: np.ndarray, base: np.ndarray, label: np.ndarray) -> np.ndarray:
    """The label of the interval that contains each depth (top inclusive, base exclusive); intervals must not overlap; a depth outside every interval gets None."""
    order = np.argsort(top)
    t, b, lab = top[order], base[order], label[order]
    idx = np.searchsorted(t, depth, side="right") - 1
    out = np.full(len(depth), None, dtype=object)
    ok = (idx >= 0) & (depth < b[np.clip(idx, 0, len(b) - 1)])
    out[ok] = lab[idx[ok]]
    return out


def log_stats_by_facies(intervals: pd.DataFrame, curves: pd.DataFrame, columns: dict[str, str], min_n: int = 5) -> pd.DataFrame:
    """For each wellbore name and genetic facies: the thickness of the interpreted intervals and, from the log samples inside them, the number of samples and the median of porosity-like columns and the geometric mean and median of
    permeability-like columns (`columns` maps a curve column to `phi` or `perm`).  A statistic from fewer than `min_n` samples is null.  Curves are matched to the intervals by the well name key."""
    rows = []
    gi = intervals[intervals["is_genetic_facies"]]
    for (w, camp), c in curves.groupby(["well_key", "curve_set"]):
        g = gi[gi["well_key"] == w]
        if g.empty:
            continue
        c = c.sort_values("md_m")
        lab = assign_interval(c["md_m"].to_numpy(float), g["top_md_m"].to_numpy(float), g["base_md_m"].to_numpy(float), g["facies_class"].to_numpy(object))
        c = c.assign(facies=lab).dropna(subset=["facies"])
        for fac, h in c.groupby("facies"):
            row = {"well_key": w, "well_name_file": g["well_name_file"].iloc[0], "curve_set": camp, "facies": fac, "samples": len(h), "thickness_m": float(g[g["facies_class"] == fac]["thickness_m"].sum())}
            for col, kind in columns.items():
                if col not in h:
                    continue
                v = h[col].astype("float64").dropna()
                if kind == "perm":
                    v = v[v > 0]
                    row[f"{col}_geomean"] = float(np.exp(np.log(v).mean())) if len(v) >= min_n else np.nan
                    row[f"{col}_median"] = float(v.median()) if len(v) >= min_n else np.nan
                else:
                    row[f"{col}_median"] = float(v.median()) if len(v) >= min_n else np.nan
            rows.append(row)
    return pd.DataFrame(rows)


def core_by_facies(plug: pd.DataFrame, intervals: pd.DataFrame, min_n: int = 5) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Each plug with a shifted depth assigned the genetic facies of its wellbore name at that depth (`facies` null outside the interpreted interval), and the plug statistics per facies: plugs, porosity median, permeability
    geometric mean and median (horizontal, Klinkenberg-corrected), blank below `min_n`.  Returns (plugs with facies, statistics).  A plug outside the interpreted interval is not assigned."""
    p = plug[plug["depth_shifted_md_m"].notna()].copy()
    p["well_key"] = p["well_name"].map(name_key)
    gi = intervals[intervals["is_genetic_facies"]]
    p["facies"] = None
    for w, g in p.groupby("well_key"):
        ig = gi[gi["well_key"] == w]
        if ig.empty:
            continue
        p.loc[g.index, "facies"] = assign_interval(g["depth_shifted_md_m"].to_numpy(float), ig["top_md_m"].to_numpy(float), ig["base_md_m"].to_numpy(float), ig["facies_class"].to_numpy(object))
    rows = []
    for (w, fac), g in p.dropna(subset=["facies"]).groupby(["well_name", "facies"]):
        por = g["por_h_pct"].astype("float64").dropna()
        k = g["kl_h_md"].astype("float64").dropna()
        k = k[k > 0]
        rows.append({"well_name": w, "facies": fac, "plugs": len(g), "porosity_median_pct": float(por.median()) if len(por) >= min_n else np.nan, "perm_geomean_md": float(np.exp(np.log(k).mean())) if len(k) >= min_n else np.nan,
                     "perm_median_md": float(k.median()) if len(k) >= min_n else np.nan, "perm_p10_md": float(k.quantile(0.1)) if len(k) >= min_n else np.nan, "perm_p90_md": float(k.quantile(0.9)) if len(k) >= min_n else np.nan})
    return p, pd.DataFrame(rows)


def dip_depth_unit_check(dips: pd.DataFrame, survey: pd.DataFrame) -> dict[str, Any]:
    """Tests the depth unit of a dip workbook against the survey of its wellbore: the well deviation of each dip is compared with the survey inclination interpolated at the stated depth taken as metres and as feet converted to metres.
    The unit under which the deviation agrees better (smaller median absolute difference) is the unit of the numbers, whatever the label says."""
    s = survey.sort_values("md_m")
    out = {}
    for name, factor in (("metres", 1.0), ("feet", 0.3048)):
        md = dips["md_stated"].to_numpy(float) * factor
        inc = np.interp(md, s["md_m"].to_numpy(float), s["inclination_deg"].to_numpy(float), left=np.nan, right=np.nan)
        diff = np.abs(inc - dips["deviation_deg"].to_numpy(float))
        out[name] = float(np.nanmedian(diff)) if np.isfinite(diff).any() else np.nan
        out[name + "_inside_survey"] = int(np.isfinite(inc).sum())
    out["unit_of_numbers"] = "metres" if (np.nan_to_num(out["metres"], nan=1e9) <= np.nan_to_num(out["feet"], nan=1e9)) else "feet"
    return out


def circular_mean_deg(angle_deg: pd.Series, weights: pd.Series | None = None) -> float:
    """The mean direction of angles in degrees (0-360), computed on the unit circle."""
    a = np.deg2rad(angle_deg.to_numpy(float))
    w = np.ones(len(a)) if weights is None else weights.to_numpy(float)
    return float(np.rad2deg(np.arctan2(np.sum(w * np.sin(a)), np.sum(w * np.cos(a)))) % 360)


def dip_bins(dips: pd.DataFrame, bin_m: float = 50.0, min_n: int = 5) -> pd.DataFrame:
    """Bed-boundary dips per depth bin (metres MD): the number of dips, the median true dip and the circular mean true dip azimuth.  Needs `md_m` (the depth in metres) and `dip_class`; bins with fewer than `min_n` dips are blank."""
    b = dips[dips["dip_class"] == "Bed Boundary"].copy()
    b["bin"] = (b["md_m"] // bin_m * bin_m).astype(float)
    rows = []
    for k, g in b.groupby("bin"):
        rows.append({"bin_top_md_m": k, "dips": len(g), "true_dip_median_deg": float(g["true_dip_deg"].median()) if len(g) >= min_n else np.nan,
                     "true_dip_azimuth_mean_deg": circular_mean_deg(g["true_dip_azimuth_deg"]) if len(g) >= min_n else np.nan})
    return pd.DataFrame(rows)
