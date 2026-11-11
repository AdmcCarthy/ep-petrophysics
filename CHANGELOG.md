# Changelog

## 1.0.0 (2026-11-11)

First public release.

## What is in it

- Checks of computed porosity, saturation, clay volume and permeability against the documented transforms.
- Rock-physics (fluid-substituted) logs compared with their inputs.
- Facies intervals and statistics, routine and special core analysis, cased-hole bond and production logs.

## Things to know

- Limits and settings that the code reads are packaged defaults and are assumptions unless their level says otherwise (see `KNOWLEDGE.md` and the defaults folder); a project file of the same name overrides them.
- Functions are pure: frames and parameters in, `(result, events)` out; nothing reads or writes files except the configuration loader.
- Licence: MIT with acknowledgement (see `LICENSE`).
- Data names from the Norwegian Offshore Directorate (Sodir, formerly NPD) appear as the sources the code was built against.
