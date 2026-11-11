# Petrophysical interpretation outputs (Volve)

Family `petrophysics_outputs` (switch in `config/pipeline_config.yaml`, section `petro_outputs` of `config/adapters/volve.yaml`, code `petrophysics`, `petrophysics`, the analysis notebooks, tests `tests/test_petrophysics.py`). Run alone with `uv run python scripts/volve_check.py --family petrophysics_outputs`.

## What is read
| Source (folder `05.PETROPHYSICAL INTERPRETATION`) | Content | Marts |
| :-- | :-- | :-- |
| `WLC_PETRO_COMPUTED_OUTPUT_1.(LAS or DLIS)` (17 wellbores) | `PHIF`, `SW`, `VSH`, `KLOGH`, `BVW` and the flags `SAND_FLAG`, `CARB_FLAG`, `COAL_FLAG` (flags from the input LAS where the output has none); 178,294 samples | `well_petro_curve` |
| `PETROPHYSICAL_REPORT_*.PDF` (text layer; 9 short reports parsed) | evaluation parameters per formation, zone results (net to gross, porosity, water saturation, permeability means), the equations and the water resistivity | `well_petro_model_param`, `_result`, `_equation`, `_fact` |
| derived | zone summaries from the logs, the arithmetic check of the permeability transform, the report zone table against the logs | `well_petro_zone_summary`, `well_petro_klogh_check`, `well_petro_zone_vs_report` |

## Rules
* The interpretation is a MODEL calibrated to core (porosity from density with a neutron correction, clay volume from the gamma ray, Archie saturation, permeability regressed on core), never a measurement.
* The report's equations are applied to the output logs as a CHECK of the arithmetic; the check assigns each sample to the first of the report's transforms that reproduces it within 5 % (assumption), a 0.001 mD floor, or none.
* The reports state no net criterion: the sand flag stands in for net (assumption); formation intervals come from the interpretation picks (the Sodir tops of the F-wells are missing); an interval that runs to the last curve sample is flagged.
* DLIS depth index `0.1 in` is converted; values are read by name; nothing is dropped.

## Findings
* One model everywhere: 43 of 51 parameter and formation combinations are identical in every report; only the gamma-ray baselines vary per well.
* The permeability log is the report's arithmetic in 7 of 8 checked wellbores (at most 1 % of samples outside the three transforms and the floor); in 15/9-F-4, 58 % of the samples (3,220-3,495 m) follow another transform (fit: log10 K = 12.1 PHIF - 9.0 VSH + 1.15, R squared 0.94) than the one the report prints (8 PHIF in the Hugin).
* Where the interval is between picks, the logs reproduce the report's mean porosity to 0.002 and water saturation to 0.06; net to gross differs by up to 0.10.
* The clay-volume equation is printed reversed in the reports (most likely typographical).

Not read: the model document the reports cite (3781-06), the longer reports (15/9-F-1, F-1 C, F-11, F-15 D, the 2005 Sleipner Ost evaluation), the first-interpretation `.LTI` files and the input curves.