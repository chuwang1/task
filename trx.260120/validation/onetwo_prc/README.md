# ONETWO PRC Kernel Validation

This directory validates `tr_onetwo_prc_kernel.f90` against PyMak's
`calculate_cyclotron_radiation_onetwo` routine.

Run from the repository root:

```bash
/home/chuw/ai-venv/bin/python validation/onetwo_prc/compare_onetwo_prc.py
```

The script:

1. Compiles `../../bpsd/bpsd_kinds.f90`, `tr_onetwo_prc_kernel.f90`, and
   `onetwo_prc_driver.f90`.
2. Runs the Fortran driver on `input_profile.csv`.
3. Loads PyMak's ONETWO function directly from `PyMak/cyclotron_radiation.py`.
4. Writes comparison CSV files and a PNG figure under `results/`.

Default validation parameters:

```text
bt_T = 6.0
minor_radius_m = 2.72
major_radius_m = 8.03
wall_reflection = 0.8
```

For validation against a complete 3-step TR run with `MDLPR=2`, see
`TR_OUTPUT_3STEP_VALIDATION.md`.
