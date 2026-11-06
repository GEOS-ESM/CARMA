# CARMA architecture

**Audience:** coding agents and developers working on the `CARMA` repository (GEOS interface layer to the Community Aerosol and Radiation Model for Atmospheres).
**Describes:** the MAPL2 line, branch `develop`.
**Verified against:** `CARMA` @ `75870b0` (2026-11-06).
**Placement:** as in the GEOSgcm fixture (`components.yaml` @ `0f2e0fb`, key `CARMA`, `local: ./src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/@GEOSchem_GridComp/@CARMA`).
Other fixtures (for example GEOSadas, GEOSldas) may place this repository
elsewhere; paths inside this repository do not change.

Parent: the Chemistry repository, see its architecture document (`../../docs/ARCHITECTURE.md` in the fixture) (section 2.3 nested repositories, section 3 the `AERO` state and the connections to CARMA). This document does not repeat that material. Terms: *ACG* is the MAPL automatic code generator (`mapl_acg`) that turns a `*_StateSpecs.rc` file into declaration and pointer include files.

## 1. Repository map

| Path | What it is |
|---|---|
| `CARMAchem_GridCompMod.F90` | MAPL gridded component `CARMAchem_GridCompMod`; public entry `SetServices` only. Also holds the registry reader (`registry_`), `add_aero`, and two methods attached to the `AERO` state (`run_aerosol_optics`, `aerosol_activation_properties`). |
| `CARMA_GridComp.F90` | Module `CARMA_GridCompMod`: types `CARMA_GridComp`, `CARMA_Registry`; `CARMA_GridCompInitialize`, `CARMA_GridCompRun` (the microphysics call per column), `CARMA_GridCompFinalize`; private `dumpGas/Element/Group` debug printers. |
| `CARMA_UtilMod.F90` | Module `CARMA_UtilMod`: host-side processes around the microphysics: `CARMA_Emissions`, `CARMA_DryDeposition`, `CARMA_WetRemoval`, `CARMA_Convection`, `CARMA_ComputeDiags`, `CARMA_GetMieTables`, `CARMA_DestroyMieTables`, `distribute_point_emissions`. |
| `ut_CARMA.F90` | A stand-alone `Program ut_CARMA` (ESMF/MAPL driver that fills import and internal states and calls `SetServices`). |
| `CARMAchem_StateSpecs.rc` | ACG specification of import, internal and export states (category blocks `IMPORT`, `INTERNAL`, `EXPORT`). |
| `CARMAchem_Registry.rc` (+ `.MIXED`, `.SU`) | CARMA configuration: bins, groups, elements, gases, process switches, optics files. Only the plain file is installed (section 5). |
| `CARMAchem_MieRegistry.rc` | Chemistry-registry style file read to build Mie tables (`Chem_RegistryCreate(rcfile='CARMAchem_MieRegistry.rc')`). |
| `CARMAchem_GridComp_ExtData.rc`, `.yaml` | ExtData entries (emissions, `CARMA_PSO4TOT`, `CARMA_HNO3`, `CARMA_H2SO4`, and similar); old-format and YAML versions are both installed. |
| `CARMA/` | Upstream CARMA library (version 3.0 per its `README`, F90 port of 2.3), embedded as a plain directory in this repository. See section 6. |
| `CMakeLists.txt` | Library build (section 6). |

No `docs/` directory existed before these files; `README.md` is a one-line title.

## 2. Purpose and place in the hierarchy

CARMA is a bin microphysics model for aerosol and cloud particles. This component runs it inside GEOS Chemistry: it keeps particle-bin and gas tracers as MAPL internal-state fields, applies emissions, microphysics, dry and wet removal and convection each Chemistry step, and exports diagnostics and an `AERO` state for radiation and moist processes (`CARMAchem_GridCompMod.F90`, SetServices, Run_).

- Parent: `GEOSchem_GridComp`, which adds it as child `CARMA` when `ENABLE_CARMA:` is true (`GEOS_ChemGridComp.F90`, SetServices), with one phase (parent document, section 2.3). The parent also makes the connections to CARMA (`AIRDENS`, `CN_PRCP`, `NCN_PRCP` from `CHEMENV`; `CARMA_PSO4TOT` from GOCART2G; `CARMA_PSOA_*` from ACHEM; `CARMA_H2SO4`, `CARMA_HNO3` from GMICHEM; `CARMA_SUSAREA`, `CARMA_SUREFF`, `CARMA_SUSAREAv` to GMICHEM or STRATCHEM when CARMA is the `AERO` provider) (`GEOS_ChemGridComp.F90`, SetServices).
- Children: none (no `MAPL_AddChild`).
- Offline/stand-alone use: `ut_CARMA.F90` (`Program ut_CARMA`) imports `SetServices` from this module.

## 3. Entry points

- **SetServices** (`CARMAchem_GridCompMod.F90`): allocates the private state `CARMAchem_State`, calls `registry_` to read `CARMAchem_Registry.rc` (and `registry_print_`), and creates a `Chem_Registry`. Only if `doing_CARMA` is true does it register Initialize, Run and Finalize with `MAPL_GridCompSetEntryPoint` (one run phase) and store the state with `ESMF_UserCompSetInternalState` under `CARMA_state`; otherwise it prints "NOT ACTIVE, defaulting to Generic No-op stubs". It then includes the ACG import spec, adds internal specs (section 4), adds export specs `AERO`, `AERO_ACI`, `AERO_DP`, includes the ACG export spec, adds timers and calls `MAPL_GenericSetServices`.
- **Initialize_**: `MAPL_GenericInitialize`; `extract_` (unpacks the private state, run-alarm interval as `cdt`, date and time); grid sizes and longitudes/latitudes in degrees; maps internal-state fields into `qa`; bootstraps the "old" fields when restart values are negative; calls `CARMA_GridCompInitialize`; reads `SULFURIC_ACID_SOURCE:` (default `tendency`); sets the `ScavengingFractionPerKm` attribute on each particle field from `FSCAV`; calls `CARMA_GetMieTables`; then builds the `AERO` state (section 4).
- **Run_** (single phase): `extract_`, then in order `CARMA_Emissions`, `CARMA_GridCompRun`, `CARMA_DryDeposition`, `CARMA_WetRemoval`, `CARMA_Convection`, `CARMA_ComputeDiags` (all take `gcCARMA, qa, import, export, nymd, nhms, cdt`).
- **Finalize_**: `CARMA_GridCompFinalize`, `CARMA_DestroyMieTables`, frees point-emission arrays and the private state (`registry_destroy_`), then `MAPL_GenericFinalize`.
- **CARMA_GridCompInitialize** (`CARMA_GridComp.F90`): creates the library object `CARMA_Create(gcCARMA%carma, NBIN, NELEM, NGROUP, ...)`, defines groups, elements and gases from the registry, adds growth (`CARMA_AddGrowth`, sulfate/H2SO4), nucleation (`CARMA_AddNucleation`, sulfate and, with a mixed group, dust/organic carbon/seasalt into the mixed group) and coagulation (`CARMA_AddCoagulation`, Fuchs kernel), then `CARMA_Initialize(...)` with the `do_*` switches from the registry. It also prepares the dust emission size fractions from the group definition.
- **CARMA_GridCompRun**: for each horizontal column (not a method on the full 3D field): gets imports and exports, optionally copies `CARMA_H2SO4`/`CARMA_HNO3` into the gas tracers (see section 5), creates a `CARMASTATE` (`CARMASTATE_Create`, or `CARMASTATE_CreateFromReference` when `do_fixedinit`), sets bins and gases (`CARMASTATE_SetBin/SetGas`), steps it (`CARMASTATE_Step`), reads bins, gases and sub-step counters back (`CARMASTATE_GetBin/GetGas/Get`) into `qa` and exports, and destroys the state (`CARMASTATE_Destroy`).

## 4. State

- **Import**: `CARMAchem_StateSpecs.rc`, `category: IMPORT`: Dynamics and Moist fields (`ZLE`, `PLE`, `T`, `Q`, `RH2`, `U`, `V`, `AIRDENS`, `CNV_MFC`, `CNV_MFD`, `CNV_QC`, `PFI_LSAN`, `PFL_LSAN`), surface fields (`USTAR`, `U10M`, `V10M`, `ZPBL`, `Z0H`, `SH`, `LWI`, `FRLAND`, `WET1`, precipitation) and source fields (`CARMA_DU_SRC`, emissions, `CARMA_PSO4TOT`, `CARMA_PSOA_*`, `CARMA_H2SO4`, `CARMA_HNO3`). `CARMA_HNO3` and `CARMA_H2SO4` are read with `notFoundOK` in `CARMA_GridComp.F90` (CARMA_GridCompRun), and `CARMA_PSO4TOT` and `CARMA_HNO3` in `CARMA_UtilMod.F90` (`CARMA_Emissions`). Declarations come from ACG includes `CARMAchem_Import___.h`, `_DeclarePointer___.h`, `_GetPointer___.h` (`CMakeLists.txt`, `mapl_acg`).
- **Internal**: not in the spec file (its `INTERNAL` block is empty); added in SetServices by `MAPL_AddInternalSpec` loops, count and names set by the registry: one tracer `CARMA::<group>::<element>::<bin>` per element and bin; one tracer per gas; then `t_old`, `<gas>_old`, `satliq_<gas>_old`, `satice_<gas>_old` (exported via `ADD2EXPORT`). Total count is `nq = NBIN*NELEM + 4*NGAS + 1` (`registry_`). Particle and gas tracers are `FRIENDLYTO='DYNAMICS:TURBULENCE'`; the "old" fields are not friendly. Initialize aborts via `_ASSERT` if the internal-spec size differs from `nq`. *Inferred:* these internal fields are checkpointed (default for `MAPL_AddInternalSpec`; no `RESTART` option is passed).
- **Export**: `AERO` (state item), `AERO_ACI` and `AERO_DP` (placeholders "not currently filled in by CARMA") in SetServices; all other exports from the ACG `EXPORT` block: per species (`DU`, `SS`, `SU`, `SM`/`OC`, `BC`, `ASH`, plus mixed-group `MX*`) fluxes, column mass, optical thickness, mass mixing ratio, surface area and effective radius, plus `CARMA_SUBSTEPS`, `CARMA_RETRIES`, `CARMA_ZSUBSTEPS`, `CARMA_SATH2SO4`. Key consumers (parent): `CARMA_SUSAREA`, `CARMA_SUREFF`, `CARMA_SUSAREAv` to GMI/StratChem.
- **`AERO` state** (built in Initialize_ only if at least one concentration element exists; `n_aerosols > 0`): a field bundle `AEROSOLS` holding references to the group concentration (`f_ienconc`) element bins in the internal state (`MAPL_FieldCreate`, `MAPL_FieldBundleAdd`); Mie table built from `CARMAchem_Registry.rc` and `CARMAchem_MieRegistry.rc` into module variable `carmaMieTable(instanceComputational)`, with attribute `mie_table_instance`; attributes `implements_aerosol_optics_method` (true), air-pressure / relative-humidity names (`PLE`, `RH`), band and wavelength 0, and field-name attributes for `EXT`, `SSA`, `ASY`; helper `add_aero` adds fields for those and for activation inputs (`T`, `FRLAND`, `SIGMA`, `NUM`, `DGN`, `density`, `KAPPA`, `FDUST`, `FSOOT`, `FORGANIC`); aerosol-cloud-interaction attributes (`aerosol_modes`, `aerosol_rmin`, `aerosol_rmrat`, `number_of_aerosol_modes`, `max_q_clean` from `MAXCLEAN:` default `1.0e-9`, `ccn_tuning` from `CCNTUNING:` default `1.8`, `cldmicro` from `CLDMICR_OPTION:` default `BACM_1M`) (`CARMAchem_GridCompMod.F90`, Initialize_).
- **Aerosol-optics methods**: `ESMF_MethodAdd(aero, label='run_aerosol_optics', userRoutine=run_aerosol_optics)` and `label='aerosol_activation_properties'` (`CARMAchem_GridCompMod.F90`, Initialize_). Executors and the wider label story: parent document (`../../docs/ARCHITECTURE.md` in the fixture) section 3.
- **Private Fortran state**: `CARMAchem_State` (`chemReg`, `gcCARMA`, `CARMAReg`, `qa`), stored under `CARMA_state`; `gcCARMA` is allocated lazily in `extract_`.

## 5. Configuration that changes structure

- Registry file name: `CARMAchem_Registry.rc` (default `rcfilen` of `CARMA_Registry`; `CARMA_GridComp.F90`); read at SetServices by `registry_`, again at Initialize for Mie tables. Keys that change structure: `NBIN`, `NGROUP`, `NELEM`, `NGAS`, `NSOLUTE`, `NWAVE` (sizes of all tracers and arrays), group and element tables (`GROUPNAME`, `IGROUP`, `ELEMNAME`, `ITYPE`, `ICOMPOSITION`, ...), `DISTRIBUTION<n>` (emission mass distributions), `GASNAME`/`IGCOMP`/`IGVAPREQ`, process switches `DO_COAG`, `DO_GROW`, `DO_SUBSTEP`, `DO_THERMO`, `DO_VDIFF`, `DO_VTRAN`, `DO_FIXEDINIT` (all optional labels), and sub-step limits (`MINSUBSTEPS`, `MAXSUBSTEPS`, `MAXRETRIES`) (`CARMAchem_GridCompMod.F90`, `registry_`).
- Group and element names matter: `CARMA_GridComp.F90` matches `DUST`, `ORGANICCARBON` and similar names in upper case, and the registry comment says "Groups/elements named 'seasalt' and 'dust' have special handling" (`CARMAchem_Registry.rc`).
- **Variants**: the shipped `CARMAchem_Registry.rc` has `NBIN 22`, 5 groups (sulfate, dust, seasalt, organiccarbon, blackcarbon), 3 gases (H2O, H2SO4, HNO3). `CARMAchem_Registry.rc.SU` has 4 groups (sulfate, dust, seasalt, smoke); `.MIXED` has `NBIN 24` and 2 groups (sulfate, `mixedp`) with 3 elements. Both variants are 483-line files that also embed old-style `ImportSpec`/`ExportSpec` blocks; neither is listed in `CMakeLists.txt` `resource_files`. *Inferred:* they are alternatives an experiment copies over `CARMAchem_Registry.rc`; no code selects among them.
- `SULFURIC_ACID_SOURCE:` (read from the component config at Initialize_, default `tendency`): `tendency` adds `CARMA_PSO4TOT*dt` to the H2SO4 gas tracer and sets the HNO3 gas tracer from the `CARMA_HNO3` import (`CARMA_UtilMod.F90`, `CARMA_Emissions`); `full_field` makes `CARMA_GridCompRun` take `CARMA_H2SO4`/`CARMA_HNO3` (mol/mol) as the gas tracers and return the updated values (`CARMA_GridComp.F90`, CARMA_GridCompRun). Who sets the key: parent document, Observation on shared config.
- ExtData: `CARMAchem_GridComp_ExtData.yaml` / `.rc` (placement of emission and chemistry sources).
- Parent flags: `ENABLE_CARMA`, `AERO_PROVIDER: CARMA` (parent document).

## 6. Dependencies and the scientific code

- Libraries (`CMakeLists.txt`): `Chem_Shared`, `Chem_Base`, `GMAO_mpeu`, `ESMF::ESMF`. From Chem_Base/Chem_Shared: `Chem_Mod`, `Chem_UtilMod`, `m_inpak90`, `Chem_Mie`, `DryDepositionMod`, `WetRemovalMod`, `ConvectionMod`, `DustEmissionMod`, `SeasaltEmissionMod`.
- Build: one library from `*.F90`/`*.F`/`*.c` globbed in `.`, `CARMA` and `CARMA/source/base`; `CARMA/source/base` is also a public include directory (module files); ACG runs on `CARMAchem_StateSpecs.rc` with `IMPORT_SPECS EXPORT_SPECS INTERNAL_SPECS GET_POINTERS DECLARE_POINTERS`; the three `.rc` files, two ExtData files and the state-spec file are installed to `etc` (`CMakeLists.txt`).
- **`CARMA/source/base`** (about 90 `.F90`/header files): the CARMA library proper, used by the interface layer through modules `carma_mod`, `carma_types_mod`, `carma_enums_mod`, `carma_constants_mod`, `carma_precision_mod`, `carmaelement_mod`, `carmagroup_mod`, `carmagas_mod`, `carmasolute_mod`, `carmastate_mod`. Files are process routines by name (for example `coagl`, `growevapl`, `sulfnuc`, `vertadv`, `step`, `mie`). Not read in detail.
- **`CARMA/tests`**, **`CARMA/doc`**, **`CARMA/bin`**: upstream stand-alone test drivers (`carma_*test.F90`, IDL readers), HTML docs and change log, a documentation tool. Not built by GEOS. The upstream `Makefile`, `make-carma.csh`, `run-*.csh` build and run the library stand-alone (`CARMA/README`); that is not the GEOS build.
- Variants of the physics are selected by registry content (what groups, gases and switches exist), not by CMake options.

## 7. Known traps

- Tracer count and order are positional: `qa(iq)` index arithmetic (`NBIN*NELEM + ...`) is shared by `CARMAchem_GridCompMod.F90`, `CARMA_GridComp.F90`, `CARMA_UtilMod.F90`; changing `NGAS` changes every offset.
- The registry is read at SetServices (before Initialize), so the registry file must be in the run directory at construction time.
- `AERO` is only filled if a concentration element exists; otherwise the exported `AERO` state has no attributes or methods.
- `AERO_ACI` and `AERO_DP` are exported but not filled.
- Optics tables: `CARMAchem_MieRegistry.rc` and the registry's `filename_optical_properties_*` entries are site paths in the shipped file (`/home/pcolarco/...`, `/discover/nobackup/...`); experiments are expected to override them (*Inferred*).
- The `SULFURIC_ACID_SOURCE` value is decided by the parent (parent document), not by this repository's files.
- Internal "old" fields use negative sentinel values to detect a missing restart (`Initialize_`).

## 8. Decisions

See [`decisions.md`](decisions.md).

## 9. Observations (noted, not fixed)

1. `CMakeLists.txt` globs `*.F90` in `.`, which includes `ut_CARMA.F90` (a `Program` unit) in the library source list; no exclusion is visible. Is the program unit meant to be built into the library?
2. In `Finalize_`, `type(Chem_Array), pointer :: qa(:)` is a local variable declared under the `!OUTPUT PARAMETERS` comment header (it is not in the argument list) and is set by `extract_` (`CARMAchem_GridCompMod.F90`, Finalize_). Neutral note.
3. `Finalize_` sets `Iam = trim(COMP_NAME) // 'Finalize_'` (no `::`), unlike the other routines.
4. `Finalize_` calls `CARMA_DestroyMieTables(gcCARMA, rc)` and does not check `rc` before continuing.
5. The `AERO` bundle comment in SetServices says it "will contain basically the same as the internal state"; the bundle created in Initialize holds only the concentration-element bins.
6. In `registry_`, `DO_FIXEDINIT` is stored from `ier(6)`, the same slot as `DO_VTRAN` (a status-array slot; noted).
7. Stale or site-specific content: shipped registry optics paths point at personal/site directories.
8. Question for maintainers: `CARMAchem_Registry.rc.MIXED` and `.SU` embed old-style spec blocks; are they still used by any experiment?
