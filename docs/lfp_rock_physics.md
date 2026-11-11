# Rock-physics (LFP) logs (Volve 15/9-19 A, BT2, SR)

Family `lfp_rock_physics` (switch in `config/pipeline_config.yaml`, section `lfp` of `config/adapters/volve.yaml`, code `lfp`, `lfp`, the analysis notebooks, tests `tests/test_lfp.py`). Run alone with `uv run python scripts/volve_check.py --family lfp_rock_physics`.

## What is read
| Source (folder `06.LFP`) | Content | Marts |
| :-- | :-- | :-- |
| `159-19A_LFP.las`, `159-19BT2_LFP.las`, `159-19SR_LFP.las` | about 170 `LFP_` curves per well (measured and edited logs, petrophysical results, velocities, impedances, Poisson's ratio and the brine / oil / gas cases `_B`, `_O`, `_G`); the selected curves of the adapter, 50,117 samples | `well_lfp_curve` |
| `lfp_petrophysics_no_15-9-19*.pdf` (interpreter's listing, 2002, text layer) | units and depth ranges, evaluation interval, shale-volume parameters, matrix density, salinity, temperature and pressure gradients (15/9-19 A), oil API, gas-oil ratio, gas gravities | `well_lfp_document`, `well_lfp_unit`, `well_lfp_shale_volume` |
| derived | physical and arithmetic checks, brine / oil / gas ordering, fluid effect in the Hugin net sand | `well_lfp_physical_check`, `well_lfp_fluid_ordering`, `well_lfp_fluid_effect` |

## Rules
* The substituted curves are a MODEL (Gassmann-type substitution with assumed fluids), not measurements; the unsubstituted elastic properties are a re-expression of the sonic and density logs.
* The data section is read as one numeric stream (records wrap over lines); null markers `-999.25` and `-9999` are set to null; implausible values are flagged and counted, never removed.
* The sand flag stands for net sand (assumption); the Hugin interval is the interpreter's listing, not the Sodir tops.
* The meaning of the suffixes `_V`, `_LOG`, `_SYNT`, `_GC` is not stated in the files.

## Findings
* Impedance = density x velocity and velocity = 304,800 / sonic slowness hold for all samples of 15/9-19 A and BT2.
* In the Hugin net sand, oil instead of brine lowers the P-impedance by 5.9-6.3 %, gas by 10.0-10.7 %; the velocity changes 3.4-4.5 %, the density supplies most of the impedance effect. Oil and gas velocities are close (the order flips in 26 % and 36 % of the samples of 15/9-19 A and BT2), the impedance order always holds.
* The interpreter assumes one oil (API 25.89, gas-oil ratio 109) for all three wellbores; it matches the 15/9-19 A laboratory oil (API 25.5, gas-oil ratio 111.8) but not 15/9-19 SR (gas-oil ratio 159.1, bubble point 273.8 bar).
* The oil saturation is set to 1 over the interval in the listing: the cases are end members.

Not read: the Word reports `no_15_9-19_a_lfp_report.doc` and `..._bt2_...doc`, the workflow and moduli, a seismic test of the substituted logs.

## Workflow definitions from the Word reports

The Word LFP evaluation reports of 15/9-19 A and 15/9-19 B (family `word_documents`, tables `word_lfp_param`, `word_lfp_top`, `word_lfp_fact`; update of January 2011 to a 1995 evaluation) define what the PDF listing left open:

- **Suffixes:** `_LOG` is the logged curve (`LFP_VP_LOG = 304800 / LFP_DT_LOG`, `LFP_AI_LOG = LFP_RHOB_LOG * LFP_VP_LOG`); `_SYNT` is a synthetic curve adjusted to the log over intervals with bad data and `must not be used for other purposes` (density from `100*(a+b*LFP_VSH)*LFP_DT_LOG**(-0.25)`, sonic by a Faust-type resistivity model); `_GC` is the shear-velocity estimate by GC; `_V` is still not defined.
- **Sand flag:** sand where porosity exceeds 0.1 and shale volume is below 0.5 (carbonate limits 77 us/ft and 2.65 g/cc); this replaces the earlier assumption that the sand flag only stands for net sand: it is the interpreter's own sand definition, still not a net-pay criterion.
- **Model flag:** `LFP_MODEL = 3 - 2 * LFP_WATER` (1 for water, 3 for hydrocarbons); in 15/9-19 A the oil interval is 3821 to 3935 m MD RKB and the evaluation runs from 3500 to 4097 m MD RKB.
- **Fluids:** the temperature and pressure equations and the oil (API 25.89, GOR 109) and gas gravities of the PDF listing are repeated and agree; salinity 150,000 ppm; matrix density 2.65 g/cc (average from core values).
- **Elevation:** the distance RKB to MSL is 78 m in both reports (water depth 85 m), against the Sodir rotary-table elevation of the wellbore: a named value for the RT/KB monitor, not applied.
- **Still open:** the fluid-substitution model itself (Gassmann-type, assumed) and the in-situ saturation behind the oil end member.

## Word report statements against the curves (mart `well_lfp_report_check`, the analysis notebooks section 5b)
Each statement of the Word LFP report (evaluation interval, oil interval, sand limits, ``LFP_SUBST = LFP_SAND``, temperature and pressure lines) is tested on the curves of the LAS file; the report does not name the porosity curve of the sand rule or the depth basis of the trend lines, so alternatives are tried and the better one is reported (``config/adapters/volve.yaml`` ``lfp.report_to_curve``, ``report_check``: tolerances are assumptions; the BT2 report file is assumed to describe the BT2 curves). Result: the oil interval (3821-3935 m MD) equals the span of model flag 3 to 0.1 m; the sand flag follows the stated limits with the effective porosity (99.8 to 100 % of samples) and less well with the total porosity; the substitution flag equals the sand flag; temperature and pressure are the stated straight lines against MD RKB (residual 0.00 to 0.09), so the pressure gradient (0.092-0.094 bar per metre MD) is not a hydrostatic gradient in a deviated well; the BT2 curve ends 20 m above the stated total depth.
