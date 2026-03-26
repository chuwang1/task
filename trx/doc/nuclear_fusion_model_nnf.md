# Nuclear Fusion Model (model_nnf) Migration

## Overview
This document records the migration of the nuclear fusion reaction switch `model_nnf` from the previous codebase to the current `trx` directory. The focus of the migration was to restore detailed fusion calculations, particularly for handling fast-ion beam-target components (NBI fast ions) which were supported by the original `TRNFDD` and `TRNFDT` subroutines but omitted in the newer `libnf` default approach.

## Key Features Migrated

1. **`model_nnf` Switch Revival:**
   - Instead of strictly using the newer `model_pnf` and `libnf` logic, the simulation can now load `model_nnf` array through standard `&TR` NAMELIST input blocks (e.g., `model_nnf(1) = 2`).
   - The array is defined in `trcomm_parm` and loaded via `trparm.f90`.
   
2. **Fast-Ion Component (`TRNFDD` & `TRNFDT`):**
   - We relocated the legacy `TRNFDD` and `TRNFDT` subroutines directly into the `trpnf.f90` module.
   - When any `model_nnf` entry is non-zero, the integration explicitly skips the default `libnf` operations and dispatches execution to either `TRNFDT` (for indices 1-4) or `TRNFDD` (for indices 11-14).
   - This reintroduces crucial fast-ion calculations like beam-target parameter (`SSB`) scaling from NBI slowing-down times (`TAUS`), generating more accurate power source tracking (`PNFCL_NSNNFNR`, `SNF_NSNNFNR`).

3. **`libsigma` Integration:**
   - The beam-target calculations in `TRNFDD/DT` depend heavily on specific empirical cross-section approximations. We introduced `libsigma.f90` to the project's build path (updated `Makefile`) to resolve the required `SIGMAM` and related helper functions.

## Implementation Details

### Modified Files:
- **`trcomm.f90`**: Reverted `model_nnf` from an allocatable array to a fixed-size integer array `INTEGER, DIMENSION(NNFM) :: model_nnf` placed in the `trcomm_parm` module so it exists before the dynamic NAMELIST loading sequence begins.
- **`trinit.f90`**: Added default initialization logic setting all elements of `model_nnf` to `0`.
- **`trparm.f90`**: Injected `model_nnf` into `NAMELIST /TR/`.
- **`trview.f90`**: Included `model_nnf(1)` to standard simulation trace logging (`WRITE(6)` outputs).
- **`trfile.f90`**: Removed a duplicate instance of `model_nnf` arrays dumping into the snapshot file which caused data misalignment on loads.
- **`libnf.f90`**: Safeguraded early-return logic `IF (ANY(model_nnf .GT. 0)) RETURN` for `model_pnf=0` instances to prevent it from resetting `nnfmax` to 0.
- **`trpnf.f90`**: Added explicit dispatch logic and embedded the missing fast-ion source routines (`TRNFDD`, `TRNFDT`, `SIGMAB`, `SIGMBS`, `SIGMAM`, `SIGMAM_DD`).
- **`Makefile`**: Included `libsigma.f90` into the `SRCS` configuration matrix.

### New Files:
- **`libsigma.f90`**: Copied from `../../../mytask/trx/` and linked via `USE libsigma`.

## Verification 
The feature was validated by loading `tr.CFEDR0114.in` where `model_nnf=2` is specified. Compilation succeeds (with `make`) and the code executes (`./tr2 < tr.CFEDR0114.in`) to successful completion (Exit Code: 0) while propagating the newly captured Fusion rates to the output CSV data dumps.
