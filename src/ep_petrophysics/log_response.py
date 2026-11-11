"""Log response per Sodir stratigraphic interval, and the stratigraphic unit of each pretest: the link between the composite logs (gamma ray, bulk density, neutron, deep resistivity, sonic, photoelectric factor) and the
wellbore formation tops (`well_strat_interval`, MD from the rig).  Pure functions on frames (`docs/borehole.md` section 10).

Rules kept here:
* tops are the official Sodir intervals (the lexicon unit id links them); the log bins are MD from the rig (RT / RKB) as the files state; no depth shift between them (a datum difference of a few metres shows up as
  a shifted boundary and is a monitored finding, not corrected);
* a statistic needs bins over at least `min_coverage` of the interval's length; below that it is blank (never a value from a few bins);
* one channel per role and wellbore, by the adapter's priority order; several files covering the same depth are reduced to their median per bin;
* the responses are measured values of the logs as delivered (composite and processed files), not interpreted properties.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

SCHEMA_VERSIONS = {"well_strat_log_response": "1.0"}

COLUMNS = {"gamma_ray": ("gr", "gAPI"), "bulk_density": ("rhob", "gcc"), "neutron": ("nphi", "vv"), "resistivity_deep": ("rdeep", "ohmm"), "sonic_slowness": ("dt", "usft"), "photoelectric": ("pef", "be")}


def one_channel_per_role(curves: pd.DataFrame, priority: dict[str, list[str]]) -> pd.DataFrame:
    """Per wellbore folder and role keep the first channel of the priority list that exists, and reduce the files that cover the same depth bin to their median.  Columns: `well_folder`, `role`, `channel`, `md_m`, `value`."""
    if curves.empty:
        return pd.DataFrame(columns=["well_folder", "role", "channel", "md_m", "value"])
    out = []
    for (folder, role), g in curves.groupby(["well_folder", "role"]):
        order = [m.upper() for m in priority.get(role, [])]
        have = {c.upper(): c for c in g["channel"].unique()}
        pick = next((have[m] for m in order if m in have), sorted(have.values())[0])
        c = g[g["channel"] == pick].groupby("md_m")["value"].median().reset_index()
        out.append(c.assign(well_folder=folder, role=role, channel=pick))
    return pd.concat(out, ignore_index=True)[["well_folder", "role", "channel", "md_m", "value"]]


def strat_response(curves1: pd.DataFrame, strat: pd.DataFrame, folder_ids: dict[str, int], min_coverage: float = 0.5) -> pd.DataFrame:
    """One row per wellbore interval of the Sodir tops (`wellbore_id`, `strat_unit`, `strat_level`, `top_md_m`, `bottom_md_m`, and the lexicon columns that exist) with, per role, the median, 10th and 90th percentile and the
    coverage (bins / interval length) from `one_channel_per_role`.  Intervals without logs keep their row with empty statistics: where there is no log is part of the answer."""
    if curves1.empty or strat.empty:
        return pd.DataFrame()
    s = strat.dropna(subset=["top_md_m", "bottom_md_m"]).copy()
    s["length_m"] = (s["bottom_md_m"].astype(float) - s["top_md_m"].astype(float)).clip(lower=1.0)
    keep = [c for c in ("wellbore_id", "strat_unit", "strat_level", "parent_unit", "unit_key", "group_name", "formation_name", "basis", "is_faulted_pick", "top_md_m", "bottom_md_m", "length_m") if c in s]
    ids = {v: k for k, v in folder_ids.items()}
    rows = []
    for wid, g in s.groupby("wellbore_id"):
        folder = ids.get(int(wid))
        c = curves1[curves1["well_folder"] == folder] if folder else curves1.iloc[0:0]
        if c.empty:
            continue
        arrays = {role: (d["md_m"].to_numpy(float), d["value"].to_numpy(float), d["channel"].iloc[0]) for role, d in c.groupby("role")}
        for r in g[keep].itertuples(index=False):
            rec = r._asdict()
            for role, (md, v, ch) in arrays.items():
                sel = (md >= rec["top_md_m"]) & (md <= rec["bottom_md_m"]) & np.isfinite(v)
                tag = COLUMNS.get(role, (role, ""))[0]
                cov = float(sel.sum()) / rec["length_m"]
                ok = cov >= min_coverage and sel.sum() >= 3
                rec[f"{tag}_median"] = float(np.median(v[sel])) if ok else np.nan
                rec[f"{tag}_p10"] = float(np.percentile(v[sel], 10)) if ok else np.nan
                rec[f"{tag}_p90"] = float(np.percentile(v[sel], 90)) if ok else np.nan
                rec[f"{tag}_coverage"] = min(cov, 1.0)
                rec[f"{tag}_channel"] = ch
            rows.append(rec)
    return pd.DataFrame(rows)


def finest_unit_at(md_m: np.ndarray, strat: pd.DataFrame) -> pd.DataFrame:
    """For each MD the finest (shortest) wellbore interval of the tops that contains it: `strat_unit`, `strat_level`, `group_name`; empty where no interval contains the depth."""
    s = strat.dropna(subset=["top_md_m", "bottom_md_m"]).copy()
    s["length_m"] = s["bottom_md_m"].astype(float) - s["top_md_m"].astype(float)
    s = s.sort_values("length_m", ascending=False)                 # longest first, so the shortest overwrites
    n = len(md_m)
    unit, level, group = np.full(n, None, dtype=object), np.full(n, None, dtype=object), np.full(n, None, dtype=object)
    for r in s.itertuples():
        hit = (md_m >= r.top_md_m) & (md_m <= r.bottom_md_m)
        unit[hit] = r.strat_unit
        level[hit] = getattr(r, "strat_level", None)
        group[hit] = getattr(r, "group_name", None)
    return pd.DataFrame({"strat_unit": unit, "strat_level": level, "group_name": group})


def attach_strat_to_pretests(pre: pd.DataFrame, strat: pd.DataFrame) -> pd.DataFrame:
    """Add the Sodir unit at each pretest depth (finest level) and, where the report names a formation, whether the report name occurs in the Sodir unit name (`report_agrees_with_sodir`; the report is the engineer's
    working name, so a disagreement is a finding to look at, not an error)."""
    if pre.empty or strat.empty:
        return pre
    out = pre.copy()
    for c in ("strat_unit", "strat_level", "group_name"):
        out[c] = None
    for wid, idx in out.dropna(subset=["wellbore_id"]).groupby("wellbore_id").groups.items():
        sw = strat[strat["wellbore_id"] == wid]
        if sw.empty:
            continue
        f = finest_unit_at(out.loc[idx, "md_m"].to_numpy(float), sw)
        for c in ("strat_unit", "strat_level", "group_name"):
            out.loc[idx, c] = f[c].to_numpy()
    out["report_agrees_with_sodir"] = [
        (str(f).lower() in str(u).lower() or str(f).lower() in str(g).lower()) if (isinstance(f, str) and isinstance(u, str)) else pd.NA for f, u, g in zip(out.get("formation", pd.Series([None] * len(out))), out["strat_unit"], out["group_name"])]
    return out


def picks_to_intervals(picks: pd.DataFrame, wid_by_name: dict[str, int], td_by_wellbore: dict[int, float] | None = None) -> pd.DataFrame:
    """Wellbore intervals from the interpretation team's picks: sorted by MD, each pick starts an interval that ends at the next pick; qualifiers NL (not logged), NR (not reached), FO (faulted out) and ER (eroded) are
    not picks; FP (faulted pick) is kept and flagged.  The result has the columns of `well_strat_interval` that exist here (`wellbore_id`, `strat_unit`, `strat_level`, `top_md_m`, `bottom_md_m`) plus `basis`.  The
    interpretation picks rank BELOW the official Sodir tops: use them only for wellbores the Sodir table does not cover."""
    p = picks[~picks["qualifier"].fillna("").isin(["NL", "NR", "FO", "ER"]) & picks["md"].notna() & (picks["surface_name"].str.lower() != "seabed")].copy()
    p["wellbore_id"] = p["well"].map(wid_by_name)
    p = p.dropna(subset=["wellbore_id"]).sort_values(["wellbore_id", "md", "obs_no"])
    rows = []
    for wid, g in p.groupby("wellbore_id"):
        md = g["md"].to_numpy(float)
        for i, r in enumerate(g.itertuples()):
            bottom = md[i + 1] if i + 1 < len(g) else (td_by_wellbore or {}).get(int(wid), np.nan)          # the last pick runs to the wellbore's TD (Sodir header), when known
            if np.isfinite(bottom) and bottom > md[i]:
                rows.append({"wellbore_id": int(wid), "strat_unit": r.surface_name, "strat_level": "interpretation pick", "top_md_m": md[i], "bottom_md_m": float(bottom), "is_faulted_pick": r.qualifier == "FP", "basis": "interpretation picks (not the official Sodir tops)"})
    return pd.DataFrame(rows)
