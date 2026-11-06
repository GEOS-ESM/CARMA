# AGENTS.md

**Audience:** coding agents and developers working on the `CARMA` repository.
**Describes:** the MAPL3 line, branch `release/MAPL-v3`.
**Verified against:** `CARMA` @ `3aa46b6` (2026-11-06).
**Placement:** as in the GEOSgcm fixture (`components.yaml` @ `f50622d`, key `CARMA`):
`src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/@GEOSchem_GridComp/@CARMA`.
Other fixtures (for example GEOSadas, GEOSldas) may place this repository
elsewhere; paths inside this repository do not change.

> **Build status on this line:** not built, not ported. On `release/MAPL-v3`, `../CMakeLists.txt` (the Chemistry `CMakeLists.txt`) builds only `GOCART` (`esma_add_subdirectories(GOCART)`, under the comment "Restore subdirs as components are ported to MAPL3."); `CARMA` appears only in a commented-out list. The Chemistry top-level library (`GEOS_ChemGridComp.F90`, which adds this component) is present but not compiled on this line: the `esma_add_library` block that builds it is commented out in `../CMakeLists.txt`. The build dependencies `Chem_Base` and `Chem_Shared` (the Chemistry `Shared` directory) are not built either. The code in this repository is the MAPL2 code, kept for a later port; it is the code of `release/MAPL-v3`, which is older than the MAPL2-line `develop` branch (`origin/develop` has 22 more commits). Differences that matter to a reader: on this line the field specs are `<ImportSpec>`/`<ExportSpec>` blocks in `CARMAchem_Registry.rc` (no `CARMAchem_StateSpecs.rc`); the group is named `smoke` (`SM` fields), not `organiccarbon` (`OC`), and black carbon is `BLACK_CARBON` in `CARMA_GridComp.F90` (`CARMA_UtilMod.F90` also tests `BLACKCARBON`); emission size distributions are fixed in `CARMA_UtilMod.F90` (no `DISTRIBUTION<n>` registry keys); the `AERO` state has no `aerosol_activation_properties` method or aerosol-cloud-interaction attributes.

## What this repository is
The GEOS interface to CARMA (Community Aerosol and Radiation Model for Atmospheres), a bin microphysics model for aerosols and clouds. One MAPL gridded component (`CARMAchem_GridCompMod`) wraps the upstream library embedded in `CARMA/`. It is a child of the Chemistry repository (`GEOSchem_GridComp`, AGENTS.md (`../AGENTS.md` in the fixture), architecture (`../docs/ARCHITECTURE.md` in the fixture); on this line the parent's documents are not written yet, the MAPL2-line versions are on its `develop` branch), which adds it when `ENABLE_CARMA` is true and makes its connections.

## Documents
| File | Content |
|---|---|
| [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) | File map, entry points, state, `AERO` state, registry configuration, build, traps, observations |
| [docs/decisions.md](docs/decisions.md) | Previous architectural decisions |
| `CHANGELOG.md` | Release notes (maintained by developers) |
| `CARMA/README`, `CARMA/doc/` | Upstream library notes (stand-alone build, not the GEOS build) |

## Where to look
- GEOS interface: `CARMAchem_GridCompMod.F90` (SetServices, Initialize_, Run_, Finalize_, registry reader, `AERO` methods); `CARMA_GridComp.F90` (library setup and per-column run); `CARMA_UtilMod.F90` (emissions, deposition, convection, diagnostics, Mie tables).
- Fields: the `<ImportSpec>` and `<ExportSpec>` blocks inside `CARMAchem_Registry.rc` (ACG input; there is no `CARMAchem_StateSpecs.rc` on this line); internal tracers are built from the registry keys in code (the `<InternalSpec>` block is empty).
- Science code: `CARMA/source/base` (the only upstream directory compiled in GEOS).

## Rules for agents
- Note bugs and oddities (document under "Observations"); do not fix them unless asked.
- Tracer indexing is positional (`NBIN*NELEM + ...`); read the "Known traps" in the architecture document before changing the registry or `NGAS`.
- Do not edit `CARMA/` (upstream) casually; GEOS-specific code is the top-level `.F90` files.
- `CARMAchem_Registry.rc.MIXED` and `.SU` are not installed by CMake; check before treating them as live.
- Keep these documents current when entry points, registry keys or connections change; update the "Verified against" line.

## Build and test
- Not built on this line (see build status above). Its `CMakeLists.txt` defines one library from sources in `.`, `CARMA`, `CARMA/source/base`, with ACG on `CARMAchem_Registry.rc`; it would be reached only through the Chemistry `CMakeLists.txt`, where it is commented out. No separate build for this repository.
- `ut_CARMA.F90` is a stand-alone driver; whether any test target runs it was not found in `CMakeLists.txt`.
- The upstream `Makefile` / `make-carma.csh` in `CARMA/` build the library stand-alone and are unrelated to GEOS builds.
