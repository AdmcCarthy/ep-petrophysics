# ep-petrophysics: What logs and cores say about the rock, and what the interpretation adds

Petrophysical interpretation outputs, rock-physics (LFP) logs, facies, log response, core analysis and special core analysis, cased-hole logs.

**Status: private, built public-ready.** The library follows the rules of a public one (no host-repository imports, no data, no private paths, an automated publication gate in `scripts/publication_gate.py`) and is published only when that gate passes. Licence: MIT with acknowledgement (`LICENSE`).

## Install

```
uv add git+https://github.com/AdmcCarthy/ep-petrophysics
```

Depends on: `ep-physics-base`, `ep-fluids`. Optional extras:

none

## What is in it

| Module | What it does |
| :--- | :--- |
| `petrophysics` | Petrophysical interpretation outputs: formation intervals, zone summaries, reproduction of the report's permeability transform and comparison with the |
| `lfp` | Rock-physics logs (LFP): physical and arithmetic checks, the effect of the pore fluid on the elastic properties, and the interpreter's assumptions aga |
| `facies` | Interpretation products: genetic facies and lithology intervals of the geomodel, the permeability logs, core by facies and the borehole-image dip inte |
| `log_response` | Log response per Sodir stratigraphic interval, and the stratigraphic unit of each pretest: the link between the composite logs (gamma ray, bulk densit |
| `core_analysis` | Conventional core analysis: cleaning, depth shift, checks and summaries.  Pure functions; they return `(result, events)` where they change or flag dat |
| `core_scal` | Special core analysis (SCAL): end points, relative permeability curves, wettability and the checks against the routine core data.  Pure functions. |
| `cased_hole` | Cased-hole logs of the Volve wells: cement bond (integrity logs) and production logging (pressure, temperature, fluid density, spinner).  Pure functio |

Functions are pure: they take frames and parameters and return `(result, events)` where they change or flag data; nothing reads or writes files except the configuration loader.

## Settings and assumptions

Every setting has a packaged default in `src/ep_petrophysics/defaults/`. A project can override a file by providing `config/<same name>` in its working directory; the project file wins. A default is an assumption unless it says otherwise: its level (`placeholder`, `inferred`, `published`, `official`, `exact`) and a one-line basis travel with it, and a citation (author and year) is given where one exists. Libraries never import a knowledge store: a host that has one passes its values in as an overlay.

## Develop

```
uv sync --all-extras
uv run pytest            # 6 test files
uv run ruff check .
uv run python scripts/publication_gate.py
```

Tests use small synthetic data generated on the fly. No network, no real data.
