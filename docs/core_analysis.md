# Core analysis (Volve 15/9-19 A and 15/9-19 BT2)

Family `core_analysis` (switch in `config/pipeline_config.yaml`, section `core` of `config/adapters/volve.yaml`, code `core_analysis`, `core_analysis`, the analysis notebooks, tests `tests/test_core_analysis.py`).
Run alone with `uv run python scripts/volve_check.py --family core_analysis`.

## What is read
| File (folder `09.CORE`) | Content | Result |
| :-- | :-- | :-- |
| `15919A.XLS`, `15919bt2.xls` | conventional (routine) core analysis, one row per plug: gas permeability, `1/Pm`, Klinkenberg-corrected permeability (horizontal, some vertical), helium porosity, retort saturations, grain density, lithology text | staging `core_plug` -> mart `well_core_plug` (860 plugs, 11 cores) |
| `* core shift.xlsx` | per core: top and base depth as cored and as shifted to the log depth (company `UNKNOWN`, method blank) | mart `well_core_shift`, shifted depth per plug |
| `05.PETROPHYSICAL INTERPRETATION/CPI/*.las` | interpretation of the same wells (`PHIF`, `KLOGH`, `SW`, `VSH`, `PORD`) | staging `core_cpi_curve`, marts `well_core_vs_log`, `well_core_vs_log_summary` |

Marts also: `well_core_summary` (per core), `well_core_klinkenberg` (implied slip factor). Not read: the core report PDFs, thin-section and rock-mechanics PDFs, the oil-water relative permeability workbooks, the `.SPW` files, `15_9-19 S` and `SR` (shift table and report only).

## Rules
* Columns are found by header name through the adapter (`core.columns`), never by position. Laboratory markers `NMP` and `NPP` are kept in `<column>_marker` and the value is null; their meaning (no measurement possible, no plug possible) is defined in the abbreviations of the final report (PDF page 13, image checked).
* No plug is removed. Plausibility flags (`flag_porosity_range`, `flag_grain_density`, `flag_saturation_sum`, `flag_klinkenberg_above_gas`, `flag_depth_order`) use assumed limits and every flagged plug stays in tables and figures (open rings).
* The raw depth is never changed; `depth_shifted_md_m` is interpolated linearly between the top and base shift of the core, `shift_basis` says `interpolated`, `end value` or `no shift table`.
* A well name is linked to a Sodir wellbore only by exact name (`link_basis`); 15/9-19 BT2 has no Sodir wellbore of that name and stays unlinked (15/9-19 B is not assumed).
* The unit of `1/Pm` is not printed; the values are consistent with the laboratory formula and with mean gas pressures of 1 to 2 bar, so 1/bar is used.
* The CPI is not independent of the core (it is normally calibrated to it): the comparison is a consistency check.

## Findings (computed in the analysis notebooks)
* **The laboratory's Klinkenberg correction is an empirical formula, not a per-plug measurement.** The final report (PDF page 11, image checked) gives it: below 2 mD `Kl = 0.68 x Kn2^1.06`, above it two iterations of `m = 0.777 x Kl^0.61`, `Kl = Kn2 - m / Pm`. Applied to the delivered Kn2 and 1/Pm it reproduces 825 of 828 corrected permeabilities (99.6 %) within 2 %; the three that differ are below 0.1 mD. (The 15/9-19 S / SR report states a correction from measurements at four back pressures, a different procedure.)
* Interpreted porosity is about 1.3-1.5 porosity units below the core (r 0.81-0.84). Interpreted permeability is 0.55-0.60 of the core median, 0.61-0.70 above the interpretation's own 0.1 mD floor (about 15 % of the interpreted values sit exactly at 0.1 mD). Causes the data cannot separate: ambient versus reservoir stress, plug versus log volume, the depth shift, calibration.
* Shifts range from -2.0 to +1.6 m; the method and company are not stated.

Primers (`config/measurement_primers.yaml`: `core_routine`, `core_depth_shift`) are from search summaries only (requests req.112, req.113).

## Measurement conditions (final reports, OCR-ed scans, key pages checked on the image)
Nitrogen permeability at steady state at **20 bar net confining sleeve pressure**; helium porosity on a Boyle's-law porosimeter and bulk volume by immersion in mercury (ambient); saturation by Dean-Stark extraction, the water saturation corrected for mud-filtrate invasion with a tritium tracer, the oil saturation not corrected (oil density 0.850 g/cc assumed); no saturation data for cores 5-7. The plug values are therefore low-stress, not reservoir-stress, values. Facts with pages are in `config/report_facts.yaml` `core_reports`.

## Special core analysis (SCAL) of 15/9-19 A
Source: the workbook `15_9_19A_o_w_relp.xls` (overview, centrifuge curves, three reservoir-condition water floods) and the 306-page report (text layer). Staging `core_scal_endpoint`, `core_scal_centrifuge`, `core_scal_waterflood`; marts `well_core_scal_endpoint`, `well_core_relperm_curve`, `well_core_scal_check_waterflood`, `well_core_scal_corey_fit`, `well_core_scal_vs_routine`. Code `core_scal`, `core_scal`, tests `tests/test_core_scal.py`.
* Three test families are kept apart (water flood at reservoir conditions on composites, centrifuge imbibition oil curves, Amott / USBM wettability on fresh, restored and cleaned plugs).
* Findings: `krw(Sor) = kw(Sor) / ko(Swi)` is above 1 in 8 of 18 tests (up to 1.48); two of three reservoir-condition curves do not start or end at the stated end points (3911 m starts at Sw 0.088, stated Swi 0.170, likely a copied starting point; 3872 m ends at Sw 0.668 against 1 - Sor 0.772); the Amott indices depend on the plug state (cleaned WWI 0.37, OWI 0; fresh WWI 0.08, OWI 0.30; restored WWI 0.22, OWI 0.15); the refitted Corey exponent is lower than the stated one in all six centrifuge plugs; SCAL permeability is about 0.8 of the routine value at the nearest plug.
* Not read: capillary pressure and resistivity index, gas displacement, gas/oil relative permeability, mercury injection (report appendices); the 15/9-19 S / SR relative permeability workbook (a garbled scan export).