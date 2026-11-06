# CARMA decisions

**Audience:** coding agents and developers working on the `CARMA` repository.
**Describes:** the MAPL3 line, branch `release/MAPL-v3`.
**Verified against:** `CARMA` @ `3aa46b6` (2026-11-06).
**Placement:** as in the GEOSgcm fixture (`components.yaml` @ `f50622d`).
Other fixtures (for example GEOSadas, GEOSldas) may place this repository
elsewhere; paths inside this repository do not change.

### Tracers declared from the registry at SetServices
- Pattern: the number and names of internal-state tracers come from `CARMAchem_Registry.rc`, read in SetServices; with no rows in the `<InternalSpec>` block (`CARMAchem_GridCompMod.F90`, SetServices, `registry_`; `CARMAchem_Registry.rc`, empty `<InternalSpec>` block).
- Likely reason: bins, groups, elements and gases are experiment choices, so the field list cannot be fixed at build time.
- Question for maintainers: none.

### Component inactive when the registry has no tracers
- Pattern: Initialize, Run and Finalize are registered only if `doing_CARMA`; otherwise MAPL generic no-op stubs remain, but exports (`AERO`, ...) are still declared (`CARMAchem_GridCompMod.F90`, SetServices).
- Likely reason: let the parent add the child unconditionally in configurations without a working registry. *Inferred.*

### Optics through the `AERO` state with Mie tables owned by this module
- Pattern: `AERO` carries attributes and fields for optics and activation, the method `run_aerosol_optics` runs `mie_` against a module-level table `carmaMieTable` (`CARMAchem_GridCompMod.F90`, Initialize_, `run_aerosol_optics`). On this line no activation-properties method is attached.
- Likely reason: radiation and moist code consume any provider's `AERO` through the same interface (see the parent's documents). *Inferred.*

### Two sources for sulfuric acid and nitric acid gas
- Pattern: `SULFURIC_ACID_SOURCE` selects `tendency` (add chemical production `CARMA_PSO4TOT*dt`) or `full_field` (use and return GMI's `CARMA_H2SO4` and `CARMA_HNO3`) (`CARMA_UtilMod.F90`, `CARMA_Emissions`; `CARMA_GridComp.F90`, CARMA_GridCompRun; `CARMAchem_GridCompMod.F90`, Initialize_).
- Likely reason: CARMA can run alone with GOCART-2G sulfate production or coupled to GMI chemistry. *Inferred.*
- Question for maintainers: the parent sets the key (see the parent's documents); why is the default `tendency`?

### Upstream library embedded as a plain directory
- Pattern: `CARMA/` is a copy of the upstream distribution (own `README`, `Makefile`, tests) inside this repository; only `CARMA/source/base` is compiled (`CMakeLists.txt`, `src_directories`).
- Likely reason: single checkout, GEOS-specific adaptation of the library. *Inferred.*
- Question for maintainers: is `CARMA/` kept in step with an upstream repository, and how?

### Per-column microphysics step
- Pattern: the library is driven one column at a time by creating, filling, stepping and destroying a `CARMASTATE` (`CARMA_GridComp.F90`, CARMA_GridCompRun).
- Likely reason: the upstream library state is a single column object. *Inferred.*
