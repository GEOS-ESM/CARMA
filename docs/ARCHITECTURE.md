# CARMA architecture

**Audience:** coding agents and developers working on the `CARMA` repository (GEOS interface layer to the Community Aerosol and Radiation Model for Atmospheres).
**Describes:** the MAPL3 line, branch `release/MAPL-v3`.
**Verified against:** `CARMA` @ `3aa46b6` (2026-11-06).
**Placement:** as in the GEOSgcm fixture (`components.yaml` @ `f50622d`, key `CARMA`, `local: ./src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/@GEOSchem_GridComp/@CARMA`).
Other fixtures (for example GEOSadas, GEOSldas) may place this repository
elsewhere; paths inside this repository do not change.

> **Build status on this line:** not built, not ported. On `release/MAPL-v3`, `../CMakeLists.txt` (the Chemistry `CMakeLists.txt`) builds only `GOCART` (`esma_add_subdirectories(GOCART)`, under the comment "Restore subdirs as components are ported to MAPL3."); `CARMA` appears only in a commented-out list (and in the commented-out `DEPENDENCIES` of the Chemistry library). The Chemistry top-level library (`GEOS_ChemGridComp.F90`, which adds this component) is present but not compiled on this line: the `esma_add_library` block that builds it is commented out in `../CMakeLists.txt`. The build dependencies `Chem_Base` and `Chem_Shared` (the Chemistry `Shared` directory, "not required by any ported component" per the parent's comment) are not built either. The code in this repository is the MAPL2 code (`MAPL_AddInternalSpec`, `MAPL_AddExportSpec`, ACG over `CARMAchem_Registry.rc`), kept for a later port. It is also older than the MAPL2-line `develop` branch (22 commits behind); statements below describe the `release/MAPL-v3` files.

Parent: the Chemistry repository, see its architecture document (`../../docs/ARCHITECTURE.md` in the fixture), which covers the nested repositories, the `AERO` state and the connections to CARMA. This document does not repeat that material. Terms: *ACG* is the MAPL automatic code generator that turns the `<ImportSpec>`/`<ExportSpec>` blocks of a resource file into declaration, spec and pointer include files.

## 1. Repository map

| Path | What it is |
|---|---|
| `CARMAchem_GridCompMod.F90` | MAPL gridded component `CARMAchem_GridCompMod`; public entry `SetServices` only. Also holds the registry reader (`registry_`), and one method attached to the `AERO` state (`run_aerosol_optics`). |
| `CARMA_GridComp.F90` | Module `CARMA_GridCompMod`: types `CARMA_GridComp`, `CARMA_Registry`; `CARMA_GridCompInitialize`, `CARMA_GridCompRun` (the microphysics call per column), `CARMA_GridCompFinalize`; private `dumpGas/Element/Group` debug printers. |
| `CARMA_UtilMod.F90` | Module `CARMA_UtilMod`: host-side processes around the microphysics: `CARMA_Emissions`, `CARMA_DryDeposition`, `CARMA_WetRemoval`, `CARMA_Convection`, `CARMA_ComputeDiags`, `CARMA_GetMieTables`, `CARMA_DestroyMieTables`, `distribute_point_emissions`. |
| `ut_CARMA.F90` | A stand-alone `Program ut_CARMA` (ESMF/MAPL driver that fills import and internal states and calls `SetServices`). |
| `CARMAchem_Registry.rc` (+ `.MIXED`, `.SU`) | CARMA configuration: bins, groups, elements, gases, process switches, optics files; also holds the ACG blocks `<ImportSpec name="CARMA">`, an empty `<InternalSpec>`, and `<ExportSpec name="CARMA">`. Only the plain file is installed (section 5). |
| `CARMAchem_MieRegistry.rc` | Chemistry-registry style file read to build Mie tables (`Chem_RegistryCreate(rcfile='CARMAchem_MieRegistry.rc')`). |
| `CARMAchem_GridComp_ExtData.rc`, `.yaml` | ExtData entries (emissions, `CARMA_PSO4TOT`, `CARMA_HNO3`, `CARMA_H2SO4`, and similar); old-format and YAML versions are both installed. The file contents differ from `develop`. |
| `CARMA/` | Upstream CARMA library (version 3.0 per its `README`, F90 port of 2.3), embedded as a plain directory in this repository. See section 6. |
| `CMakeLists.txt` | Library build (section 6); not reached on this line. |

No `docs/` directory existed before these files; `README.md` is a one-line title.

## 2. Purpose and place in the hierarchy

CARMA is a bin microphysics model for aerosol and cloud particles. This component runs it inside GEOS Chemistry: it keeps particle-bin and gas tracers as MAPL internal-state fields, applies emissions, microphysics, dry and wet removal and convection each Chemistry step, and exports diagnostics and an `AERO` state for radiation and moist processes (`CARMAchem_GridCompMod.F90`, SetServices, Run_).

- Parent: `GEOSchem_GridComp`, which adds it as child `CARMA` when `ENABLE_CARMA:` is true (`GEOS_ChemGridComp.F90`, SetServices), with one phase. The parent also makes the connections to CARMA (`AIRDENS`, `CN_PRCP`, `NCN_PRCP` from `CHEMENV`; `CARMA_PSO4TOT` from GOCART2G; `CARMA_PSOA_*` from ACHEM; `CARMA_H2SO4`, `CARMA_HNO3` from GMICHEM; `CARMA_SUSAREA`, `CARMA_SUREFF`, `CARMA_SUSAREAv` to GMICHEM or STRATCHEM when CARMA is the `AERO` provider) (`GEOS_ChemGridComp.F90`, SetServices).
- Children: none (no `MAPL_AddChild`).
- Offline/stand-alone use: `ut_CARMA.F90` (`Program ut_CARMA`) imports `SetServices` from this module.

## 3. Entry points

- **SetServices** (`CARMAchem_GridCompMod.F90`): allocates the private state `CARMAchem_State`, calls `registry_` to read `CARMAchem_Registry.rc` (and `registry_print_`), and creates a `Chem_Registry`. Only if `doing_CARMA` is true does it register Initialize, Run and Finalize with `MAPL_GridCompSetEntryPoint` (one run phase) and store the state with `ESMF_UserCompSetInternalState` under `CARMA_state`; otherwise it prints "NOT ACTIVE, defaulting to Generic No-op stubs". It then includes the ACG import spec (`CARMA_ImportSpec___.h`), adds internal specs (section 4), adds export specs `AERO`, `AERO_ACI`, `AERO_DP`, includes the ACG export spec (`CARMA_ExportSpec___.h`), adds timers and calls `MAPL_GenericSetServices`.
- **Initialize_**: `MAPL_GenericInitialize`; `extract_` (unpacks the private state, run-alarm interval as `cdt`, date and time); grid sizes and longitudes/latitudes in degrees; maps internal-state fields into `qa`; bootstraps the "old" fields when restart values are negative; calls `CARMA_GridCompInitialize`; reads `SULFURIC_ACID_SOURCE:` (default `tendency`); sets the `ScavengingFractionPerKm` attribute on each particle field from `FSCAV`; calls `CARMA_GetMieTables`; then builds the `AERO` state (section 4). The `AERO` state methods are only `run_aerosol_optics`; there is no activation-properties method on this line.
- **Run_** (single phase): `extract_`, then in order `CARMA_Emissions`, `CARMA_GridCompRun`, `CARMA_DryDeposition`, `CARMA_WetRemoval`, `CARMA_Convection`, `CARMA_ComputeDiags` (all take `gcCARMA, qa, import, export, nymd, nhms, cdt`).
- **Finalize_**: `CARMA_GridCompFinalize`, `CARMA_DestroyMieTables`, frees point-emission arrays and the private state (`registry_destroy_`), then `MAPL_GenericFinalize`.
- **CARMA_GridCompInitialize** (`CARMA_GridComp.F90`): creates the library object `CARMA_Create(gcCARMA%carma, NBIN, NELEM, NGROUP, ...)`, defines groups, elements and gases from the registry, adds growth (`CARMA_AddGrowth`, sulfate/H2SO4), nucleation (`CARMA_AddNucleation`, sulfate and, with a mixed group, dust, smoke and seasalt into the mixed group: `ielm_dust`, `ielm_smoke`, `ielm_seasalt`) and coagulation (`CARMA_AddCoagulation`, Fuchs kernel), then `CARMA_Initialize(...)` with the `do_*` switches from the registry. It also prepares the dust emission size fractions from the group definition.
- **CARMA_GridCompRun**: for each horizontal column (not a method on the full 3D field): gets imports and exports, optionally copies `CARMA_H2SO4`/`CARMA_HNO3` into the gas tracers (see section 5), creates a `CARMASTATE` (`CARMASTATE_Create`, or `CARMASTATE_CreateFromReference` when `do_fixedinit`), sets bins and gases (`CARMASTATE_SetBin/SetGas`), steps it (`CARMASTATE_Step`), reads bins, gases and sub-step counters back (`CARMASTATE_GetBin/GetGas/Get`) into `qa` and exports, and destroys the state (`CARMASTATE_Destroy`).

## 4. State

- **Import**: block `<ImportSpec name="CARMA">` in `CARMAchem_Registry.rc`: Dynamics and Moist fields (`ZLE`, `PLE`, `T`, `Q`, `RH2`, `U`, `V`, `AIRDENS`, `CNV_MFC`, `CNV_MFD`, `CNV_QC`, `PFI_LSAN`, `PFL_LSAN`), surface fields (`USTAR`, `U10M`, `V10M`, `ZPBL`, `Z0H`, `SH`, `LWI`, `FRLAND`, `WET1`, precipitation) and source fields (`CARMA_DU_SRC`, emissions, `CARMA_PSO4TOT`, `CARMA_PSOA_*`, `CARMA_H2SO4`, `CARMA_HNO3`). `CARMA_HNO3` and `CARMA_H2SO4` are read with `notFoundOK` in `CARMA_GridComp.F90` (CARMA_GridCompRun), and `CARMA_PSO4TOT` and `CARMA_HNO3` in `CARMA_UtilMod.F90` (`CARMA_Emissions`). The emission source imports are `CARMA_SM_BIOMASS`, `CARMA_SM_BIOFUEL`, `CARMA_SM_ANTEOC1`, `CARMA_SM_ANTEOC2`, `CARMA_SM_SHIP`, `CARMA_OC_TERPENE`, plus `CARMA_PSOA_ANTHRO_VOC` and `CARMA_PSOA_BIOB_VOC` (`CARMAchem_Registry.rc` `<ImportSpec>`; read in `CARMA_UtilMod.F90`, `CARMA_Emissions`). Declarations come from the ACG includes `CARMA_ImportSpec___.h`, `CARMA_ExportSpec___.h`, `CARMA_DeclarePointer___.h`, `CARMA_GetPointer___.h` (`CARMAchem_GridCompMod.F90`; `CMakeLists.txt`, `new_esma_generate_automatic_code` on `CARMAchem_Registry.rc`, which lists `CARMA_ExportSpec___.h` and `CARMA_GetPointer___.h`).
- **Internal**: not in the registry file (its `<InternalSpec>` block is empty); added in SetServices by `MAPL_AddInternalSpec` loops, count and names set by the registry: one tracer `CARMA::<group>::<element>::<bin>` per element and bin; one tracer per gas; then `t_old`, `<gas>_old`, `satliq_<gas>_old`, `satice_<gas>_old` (exported via `ADD2EXPORT`). Total count is `nq = NBIN*NELEM + 4*NGAS + 1` (`registry_`). Particle and gas tracers are `FRIENDLYTO='DYNAMICS:TURBULENCE'`; the "old" fields are not friendly. Initialize aborts via `_ASSERT` if the internal-spec size differs from `nq`. *Inferred:* these internal fields are checkpointed (default for `MAPL_AddInternalSpec`; no `RESTART` option is passed).
- **Export**: `AERO` (state item), `AERO_ACI` and `AERO_DP` (placeholders "not currently filled in by CARMA") in SetServices; all other exports from the `<ExportSpec>` block of `CARMAchem_Registry.rc`: per species (`DU`, `SS`, `SU`, `SM`, `BC`, `ASH`, plus mixed-group `MX*`; the smoke fields are `CARMA_SM*`, for example `CARMA_SMSD`, `CARMA_SMMASS`) fluxes, column mass, optical thickness, mass mixing ratio, surface area and effective radius, plus `CARMA_SUBSTEPS`, `CARMA_RETRIES`, `CARMA_ZSUBSTEPS`, `CARMA_SATH2SO4`. Key consumers (parent): `CARMA_SUSAREA`, `CARMA_SUREFF`, `CARMA_SUSAREAv` to GMI/StratChem.
- **`AERO` state** (the bundle and the attribute `implements_aerosol_optics_method` (true) are set in Initialize_ unconditionally; the rest only if at least one concentration element exists, `n_aerosols > 0`): a field bundle `AEROSOLS` holding references to the group concentration (`f_ienconc`) element bins in the internal state (`MAPL_FieldCreate`, `MAPL_FieldBundleAdd`); Mie table built from `CARMAchem_Registry.rc` and `CARMAchem_MieRegistry.rc` into module variable `carmaMieTable(instanceComputational)`, with attribute `mie_table_instance`; attributes for air-pressure / relative-humidity names (`PLE`, `RH`), band 0, a cloud-fraction name left empty, and field-name attributes for `EXT`, `SSA`, `ASY`; the fields for `PLE`, `RH`, `EXT`, `SSA`, `ASY` are created inline with `MAPL_FieldCreateEmpty` / `MAPL_FieldAllocCommit` and added with `MAPL_StateAdd` (`CARMAchem_GridCompMod.F90`, Initialize_). On this line there is no helper `add_aero`, no wavelength attribute, no activation inputs and no aerosol-cloud-interaction attributes (`aerosol_modes`, `max_q_clean`, `ccn_tuning`, `cldmicro`); `AERO_ACI` is not filled (see traps).
- **Aerosol-optics method**: `ESMF_MethodAdd(aero, label='run_aerosol_optics', userRoutine=run_aerosol_optics)` (`CARMAchem_GridCompMod.F90`, Initialize_); no other method is attached on this line. Executors and the wider label story: the parent's documents.
- **Private Fortran state**: `CARMAchem_State` (`chemReg`, `gcCARMA`, `CARMAReg`, `qa`), stored under `CARMA_state`; `gcCARMA` is allocated lazily in `extract_`.

## 5. Configuration that changes structure

- Registry file name: `CARMAchem_Registry.rc` (default `rcfilen` of `CARMA_Registry`; `CARMA_GridComp.F90`); read at SetServices by `registry_`, again at Initialize for Mie tables. Keys that change structure: `NBIN`, `NGROUP`, `NELEM`, `NGAS`, `NSOLUTE`, `NWAVE` (sizes of all tracers and arrays), group and element tables (`GROUPNAME`, `IGROUP`, `ELEMNAME`, `ITYPE`, `ICOMPOSITION`, ...), `GASNAME`/`IGCOMP`/`IGVAPREQ`, process switches `DO_COAG`, `DO_GROW`, `DO_SUBSTEP`, `DO_THERMO`, `DO_VDIFF`, `DO_VTRAN`, `DO_FIXEDINIT` (all optional labels), and sub-step limits (`MINSUBSTEPS`, `MAXSUBSTEPS`, `MAXRETRIES`) (`CARMAchem_GridCompMod.F90`, `registry_`).
- There is no `DISTRIBUTION<n>` key on this line: emission size distributions are fixed in `CARMA_UtilMod.F90` (`dMbc`, `dMpin`).
- Group and element names matter: `CARMA_GridComp.F90` matches `DUST`, `SMOKE`, `BLACK_CARBON` and similar names in upper case,. `CARMA_UtilMod.F90` tests `'blackcarbon'`/`'BLACKCARBON'` (and mixed-group element `BLACKCARBON`) while `CARMA_GridComp.F90` matches `'BLACK_CARBON'`: two spellings of the black-carbon group name in the same component (neutral note).
- **Variants**: the shipped `CARMAchem_Registry.rc` has `NBIN 22`, `NGROUP 1`, `NELEM 1`, `NGAS 0`, while its group and element tables list three entries (`GROUPNAME: dust seasalt smoke`, `ELEMNAME: pc pc pc`); the file's comment says tables "must be at least NGROUP entries", so *inferred:* only the first `NGROUP`/`NELEM` entries are used. `CARMAchem_Registry.rc.SU` has 4 groups (sulfate, dust, seasalt, smoke); `.MIXED` has `NBIN 24` and 2 groups (sulfate, `mixedp`) with 3 elements. Both variants are 483-line files that also embed old-style `ImportSpec`/`ExportSpec` blocks; neither is listed in `CMakeLists.txt` `resource_files`. *Inferred:* they are alternatives an experiment copies over `CARMAchem_Registry.rc`; no code selects among them.
- `SULFURIC_ACID_SOURCE:` (read from the component config at Initialize_, default `tendency`): `tendency` adds `CARMA_PSO4TOT*dt` to the H2SO4 gas tracer and sets the HNO3 gas tracer from the `CARMA_HNO3` import (`CARMA_UtilMod.F90`, `CARMA_Emissions`); `full_field` makes `CARMA_GridCompRun` take `CARMA_H2SO4`/`CARMA_HNO3` (mol/mol) as the gas tracers and return the updated values (`CARMA_GridComp.F90`, CARMA_GridCompRun). Who sets the key: the parent's documents.
- ExtData: `CARMAchem_GridComp_ExtData.yaml` / `.rc` (placement of emission and chemistry sources).
- Parent flags: `ENABLE_CARMA`, `AERO_PROVIDER: CARMA` (parent).

## 6. Dependencies and the scientific code

- Libraries (`CMakeLists.txt`): `Chem_Shared`, `Chem_Base` (neither built on this line), `GMAO_mpeu`, `ESMF::ESMF`. From Chem_Base/Chem_Shared: `Chem_Mod`, `Chem_UtilMod`, `m_inpak90`, `Chem_Mie`, `DryDepositionMod`, `WetRemovalMod`, `ConvectionMod`, `DustEmissionMod`, `SeasaltEmissionMod`.
- Build (not built on this line; the Chemistry `CMakeLists.txt` does not add `CARMA`): one library from `*.F90`/`*.F`/`*.c` globbed in `.`, `CARMA` and `CARMA/source/base`; `CARMA/source/base` is also a public include directory (module files); `new_esma_generate_automatic_code` runs on `CARMAchem_Registry.rc` for `CARMA_ExportSpec___.h;CARMA_GetPointer___.h`; `CARMAchem_Registry.rc`, `CARMAchem_MieRegistry.rc` and the two ExtData files are installed to `etc` (`CMakeLists.txt`).
- **`CARMA/source/base`** (about 90 `.F90`/header files): the CARMA library proper, used by the interface layer through modules `carma_mod`, `carma_types_mod`, `carma_enums_mod`, `carma_constants_mod`, `carma_precision_mod`, `carmaelement_mod`, `carmagroup_mod`, `carmagas_mod`, `carmasolute_mod`, `carmastate_mod`. Files are process routines by name (for example `coagl`, `growevapl`, `sulfnuc`, `vertadv`, `step`, `mie`). Not read in detail.
- **`CARMA/tests`**, **`CARMA/doc`**, **`CARMA/bin`**: upstream stand-alone test drivers (`carma_*test.F90`, IDL readers), HTML docs and change log, a documentation tool. Not built by GEOS. The upstream `Makefile`, `make-carma.csh`, `run-*.csh` build and run the library stand-alone (`CARMA/README`); that is not the GEOS build.
- Variants of the physics are selected by registry content (what groups, gases and switches exist), not by CMake options.

## 7. Known traps

- Tracer count and order are positional: `qa(iq)` index arithmetic (`NBIN*NELEM + ...`) is shared by `CARMAchem_GridCompMod.F90`, `CARMA_GridComp.F90`, `CARMA_UtilMod.F90`; changing `NGAS` changes every offset.
- The registry is read at SetServices (before Initialize), so the registry file must be in the run directory at construction time.
- The `AEROSOLS` bundle and `implements_aerosol_optics_method` (true) are set unconditionally in Initialize_; the Mie table, the other attributes and `ESMF_MethodAdd` (`run_aerosol_optics`) happen only if a concentration element exists (`n_aerosols > 0`).
- `AERO_ACI` and `AERO_DP` are exported but not filled (comments in SetServices say so).
- Optics tables: `CARMAchem_MieRegistry.rc` and the registry's `filename_optical_properties_*` entries (keys `DU`, `SS`, `BC`, `SM`, `SU` on this line) are site paths in the shipped file (`/home/pcolarco/...`, `/discover/nobackup/...`); experiments are expected to override them (*Inferred*).
- The `SULFURIC_ACID_SOURCE` value is decided by the parent (see the parent's documents), not by this repository's files.
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
8. In the shipped `CARMAchem_Registry.rc`, `filename_optical_properties_BC` names a `carma_optics_SU` file (same file as the `SU` entry).
9. Question for maintainers: `CARMAchem_Registry.rc.MIXED` and `.SU` embed old-style spec blocks; are they still used by any experiment?
