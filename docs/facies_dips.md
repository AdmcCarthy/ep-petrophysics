# Facies, permeability logs and image dips (Volve)

Family `facies_dips` (switch in `config/pipeline_config.yaml`, section `facies_dips` of `config/adapters/volve.yaml`, code `facies`, `facies`, the analysis notebooks, tests `tests/test_facies.py`). Run alone with `uv run python scripts/volve_check.py --family facies_dips`.

## What is read
| Source | Content | Marts |
| :-- | :-- | :-- |
| `*Facies.xlsx` (05.PETROPHYSICAL INTERPRETATION, 11 wellbore names) | geomodel facies intervals: `GENETIC FACIES` (source `GM09`, 139 intervals of the Hugin interval, 15 classes) and `LITH` (source `SFINX`, 7,856 fine lithology-class intervals of 15/9-19 A and BT2) | `well_facies_interval` |
| `*KLOGH_NEW.las` (geomodel permeability, 8 wellbores) and `*_CPI.las` (15/9-19 A, BT2, SR: PHIF, KLOGH, SW, VSH) | log curves inside the facies intervals | `well_facies_curve`, `well_facies_log_stats` |
| `Statoil_Volve_15_9-F-14_Dips.xls` (07.IMAGE) | 469 dips of an azimuthal resistivity image: bed boundaries (428), possible faults (27), unconformable surfaces (14) with deviation, hole azimuth, true and apparent dip | `well_image_dip`, `well_image_dip_bin` |
| plugs of the analysis notebooks | assigned a genetic facies at their shifted depth | `well_core_plug_facies`, `well_core_by_facies` |

## Rules
* A facies is the geomodellers' interpretation (definitions not read); intervals are shown as delivered, gaps and overlaps flagged; statistics need 5 samples; plugs are assigned only through the shifted depth and only inside an interpreted interval.
* The facies colours are a neutral categorical palette (no standard for facies; no geological colour is invented).
* The depth unit of the dip workbook is TESTED against the survey inclination (the sheet labels its depths `ft`; the numbers are metres: the deviations agree to 0.002 degrees as metres and 64 degrees as feet).
* Dip azimuths are averaged as angles; no structural-dip removal, fault throw or sealing statement is made.

## Findings
* Inside the GM09 interval, tidal channel plugs average 1,058-1,153 mD (geometric mean) against 47-102 mD for tidal bar plugs at a similar median porosity; facies alone explain 28 % of the variance of log10 permeability of 410 plugs, porosity alone 53 %, both 72 %.
* The F-14 dips: bed boundaries have a median true dip of 23 degrees (10th-90th percentile 7-42) in a well deviated 40-81 degrees; 6 'unconformable surface' picks lie within 30 m of the biostratigraphic unconformable boundaries of the same well (2,926-3,000 m MD).

Not read: facies definitions and the 2009 geomodel report, the `SFINX` classes, the image files (`.pds`, `.dlis`), the F-14 image report and presentations, the facies of 15/9-19 SR beyond 6 intervals.