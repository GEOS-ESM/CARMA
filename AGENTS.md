# AGENTS.md

**Audience:** coding agents and developers working on the `CARMA` repository.
**Describes:** the MAPL2 line, branch `develop`.
**Verified against:** `CARMA` @ `75870b0` (2026-11-06).
**Placement:** as in the GEOSgcm fixture (`components.yaml` @ `0f2e0fb`, key `CARMA`):
`src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/@GEOSchem_GridComp/@CARMA`.
Other fixtures (for example GEOSadas, GEOSldas) may place this repository
elsewhere; paths inside this repository do not change.

## What this repository is
The GEOS interface to CARMA (Community Aerosol and Radiation Model for Atmospheres), a bin microphysics model for aerosols and clouds. One MAPL gridded component (`CARMAchem_GridCompMod`) wraps the upstream library embedded in `CARMA/`. It is a child of the Chemistry repository (`GEOSchem_GridComp`, [AGENTS.md](../AGENTS.md), [architecture](../docs/ARCHITECTURE.md)), which adds it when `ENABLE_CARMA` is true and makes its connections.

## Documents
| File | Content |
|---|---|
| [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) | File map, entry points, state, `AERO` state, registry configuration, build, traps, observations |
| [docs/decisions.md](docs/decisions.md) | Previous architectural decisions |
| `CHANGELOG.md` | Release notes (maintained by developers) |
| `CARMA/README`, `CARMA/doc/` | Upstream library notes (stand-alone build, not the GEOS build) |

## Where to look
- GEOS interface: `CARMAchem_GridCompMod.F90` (SetServices, Initialize_, Run_, Finalize_, registry reader, `AERO` methods); `CARMA_GridComp.F90` (library setup and per-column run); `CARMA_UtilMod.F90` (emissions, deposition, convection, diagnostics, Mie tables).
- Fields: `CARMAchem_StateSpecs.rc` (ACG); internal tracers are built from `CARMAchem_Registry.rc` in code, not in the spec file.
- Science code: `CARMA/source/base` (the only upstream directory compiled in GEOS).

## Rules for agents
- Note bugs and oddities (document under "Observations"); do not fix them unless asked.
- Tracer indexing is positional (`NBIN*NELEM + ...`); read the "Known traps" in the architecture document before changing the registry or `NGAS`.
- Do not edit `CARMA/` (upstream) casually; GEOS-specific code is the top-level `.F90` files.
- `CARMAchem_Registry.rc.MIXED` and `.SU` are not installed by CMake; check before treating them as live.
- Keep these documents current when entry points, registry keys or connections change; update the "Verified against" line.

## Build and test
- Built as a library by the GEOSgcm CMake build (`CMakeLists.txt`): sources from `.`, `CARMA`, `CARMA/source/base`; ACG on `CARMAchem_StateSpecs.rc`. No separate build for this repository.
- `ut_CARMA.F90` is a stand-alone driver; whether any test target runs it was not found in `CMakeLists.txt`.
- The upstream `Makefile` / `make-carma.csh` in `CARMA/` build the library stand-alone and are unrelated to GEOS builds.
