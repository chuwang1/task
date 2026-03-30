# MODELG5 Non-Symmetric LCFS Axis Issue

## Background

We investigated the `MODELG=5` path in `TASK/EQ` to understand whether a non-symmetric gfile LCFS is preserved and why the final magnetic axis still drifts toward low `Z` even after the LCFS shape becomes visibly up-down asymmetric.

The original concern was:

- the gfile LCFS is non-symmetric;
- `TASK/EQ` appears to read it;
- but later outputs often look closer to a symmetric or low-`Z` solution.

## Confirmed Code Path

### 1. Original gfile LCFS is read

The gfile reader loads LCFS points into:

- `RSU(NSUM+1)`
- `ZSU(NSUM+1)`

References:

- `eq/eqcom1.inc:47`
- `eq/eq-eqdsk.f:49`
- `eq/eq-eqdsk.f:50`

This means the original boundary points from the gfile do enter the code.

### 2. MODELG=5 boundary setup uses gfile LCFS

`EQDEFB` uses the gfile LCFS path when:

- `MODELG=5`
- `NSUMAX>=8`
- `IUSELCFS!=0`

References:

- `eq/eqcalc.f:216`
- `eq/eqcalc.f:217`

The conversion routine is:

- `EQSET_RHOB_FROM_RSU`

It converts `RSU/ZSU` into internal boundary arrays:

- `RHOM`
- `RHOG`

References:

- `eq/eqcalc.f:311`
- `eq/eqcalc.f:326`
- `eq/eqcalc.f:327`
- `eq/eqcalc.f:399`
- `eq/eqcalc.f:404`

### 3. We fixed the main LCFS-loss point in EQTORZ

Originally, `EQTORZ` rebuilt the boundary radius from analytic `EQFBND`, which depends only on symmetric parameterized geometry (`RA`, `RKAP`, `RDLT`) and therefore can destroy true gfile asymmetry.

Original logic reference:

- `eq/eqcalc.f:1077` to `eq/eqcalc.f:1091`

We changed this so that for `MODELG=5` with `IUSELCFS!=0`, `EQTORZ` now uses interpolated `RHOG` first, and only falls back to `EQFBND` if interpolation fails.

Related additions:

- `EQGET_RHOB` in `eq/eqcalc.f:441`
- `EQTORZ` branch in `eq/eqcalc.f:1077`

## What Is Working Now

After the `EQTORZ` change:

- the exported / reconstructed LCFS shape becomes visibly up-down asymmetric;
- this confirms the previous main geometry-loss mechanism was indeed `EQFBND` inside `EQTORZ`.

So the non-symmetric boundary is now preserved better than before.

## Current Problem

Even though the LCFS shape is now asymmetric, the final magnetic axis still drifts toward lower `Z`, often ending up near `Z=0` or below the gfile axis location.

In other words:

- boundary asymmetry is preserved better;
- axis position is still not preserved.

This suggests the remaining problem is no longer mainly the LCFS representation, but the way the equilibrium iteration and axis search behave after `PSIRZ` is rebuilt.

## What We Tried

### Attempt A: Shift internal grid center from `(RR,0)` to `(RAXIS,ZAXIS)`

We tested changing the internal mesh and `EQTORZ` polar reference from:

- center at `RR,0`

to:

- center at `RAXIS,ZAXIS`

Modified lines included:

- `eq/eqcalc.f:262`
- `eq/eqcalc.f:268`
- `eq/eqcalc.f:274`
- `eq/eqcalc.f:275`
- `eq/eqcalc.f:1075`
- `eq/eqcalc.f:1091`

Observed result:

- the magnetic axis did move away from `Z=0`;
- but the solver became unstable;
- many `XX EQMAGS: NOT ENOUGH N` and `XX EQAXIS: NEWTN ERROR` appeared;
- final output file was not reliably produced.

Conclusion:

- this proves the axis position is sensitive to the internal reference center;
- however, a direct center shift is not a safe minimal fix in the current code structure.

This change was reverted.

### Attempt B: Change EQMAGS closure criterion

We changed `EQMAGS` so that it no longer detects one closed turn by checking whether the orbit crosses `Z=ZINIT` again.

Original logic references:

- `eq/eqsub.f:161`
- `eq/eqsub.f:163`
- `eq/eqsub.f:191`
- `eq/eqsub.f:201`

New logic:

- accumulate poloidal angle around `(RAXIS,ZAXIS)`;
- require accumulated angle about `2*pi`;
- require return near the starting point.

Modified region:

- `eq/eqsub.f:158` to `eq/eqsub.f:204`

Observed result:

- code still compiles and runs;
- output file is produced;
- but `XX EQMAGS: NOT ENOUGH N` still appears many times;
- axis drift behavior is not significantly improved.

Conclusion:

- the old horizontal-line closure condition was not ideal;
- but it is not the main reason the axis is pulled downward.

### Attempt C: Change EQAXIS bounds check to use LCFS bounding box

We modified `EQAXIS` so that, for `MODELG=5` and `IUSELCFS!=0`, axis bounds are checked against the real `RSU/ZSU` bounding box rather than the old symmetric box based on `RB` and `RKAP`.

Modified region:

- `eq/eqsub.f:58` to `eq/eqsub.f:81`

Observed result:

- code still compiles and runs;
- output file is produced;
- axis behavior remains almost unchanged.

Conclusion:

- the symmetric box in `EQAXIS` is not the main cause either.

## Current Best Understanding

At this stage, the main findings are:

1. The original gfile LCFS is read correctly.
2. `EQDEFB` can use that LCFS and convert it to `RHOM/RHOG`.
3. The previous major asymmetry-loss point in `EQTORZ` has been fixed by using `RHOG` instead of `EQFBND`.
4. The LCFS shape now looks non-symmetric as expected.
5. The remaining unresolved issue is the magnetic axis location, not the boundary shape.

The most likely remaining problem is:

- the equilibrium field `PSIRZ` produced during the iteration develops a lower-`Z` extremum structure;
- then `find_axis` / Newton follows that structure consistently.

This means the root cause is probably in one of these areas:

- how `PSI(NTG,NSG)` evolves inside `EQLOOP`;
- how `PSI(sigma,theta)` is mapped back into `PSIRZ` by `EQTORZ`;
- how `PSIGD` and `find_axis` see the local gradient structure after remapping.

## Most Relevant Current Files and Locations

### LCFS reading and storage

- `eq/eqcom1.inc:47`
- `eq/eq-eqdsk.f:49`
- `eq/eq-eqdsk.f:50`

### MODELG=5 LCFS boundary conversion

- `eq/eqcalc.f:216`
- `eq/eqcalc.f:217`
- `eq/eqcalc.f:311`
- `eq/eqcalc.f:399`
- `eq/eqcalc.f:404`

### Stable non-symmetric boundary preservation in EQTORZ

- `eq/eqcalc.f:441`
- `eq/eqcalc.f:1077`
- `eq/eqcalc.f:1078`
- `eq/eqcalc.f:1085`

### Axis search

- `eq/eqsub.f:5`
- `eq/eqsub.f:26`
- `eq/eqsub.f:41`
- `eq/eqsub.f:91`
- `eq/eqsub.f:270`
- `eq/eqsub.f:284`

### Field-line / magnetic-surface tracing

- `eq/eqsub.f:116`
- `eq/eqsub.f:138`
- `eq/eqsub.f:158`
- `eq/eqsub.f:192`

### PSIRZ-based radial surface tracing in EQCALQ

- `eq/eqcalq.f:194`
- `eq/eqcalq.f:195`
- `eq/eqcalq.f:203`

## Recommended Next Step

The next debugging step should not focus on LCFS geometry anymore.

Instead, it should directly diagnose why the axis of the rebuilt `PSIRZ` moves downward.

Recommended diagnostics:

1. In `EQAXIS`, print `PSI` or `grad PSI` above and below the current axis guess.
2. After `EQTORZ`, scan `PSIRZ` for the actual minimum / axis candidate location on the `R,Z` grid.
3. Compare that location with the Newton result from `find_axis`.

This will tell us whether the downward drift comes from:

- the field itself (`PSIRZ` already has the wrong axis), or
- the axis finder (`find_axis` / `PSIGD` / Newton chooses the wrong extremum).

## Current Status Summary

- stable improvement kept: `EQTORZ` now preserves non-symmetric LCFS using `RHOG` interpolation;
- unstable attempt reverted: shifting internal grid center to `RAXIS/ZAXIS`;
- limited-gain attempts kept for testing so far:
  - `EQMAGS` new closure criterion;
  - `EQAXIS` LCFS-based bounding box.

The project is now at the stage where the main unresolved issue is:

> Why does the rebuilt equilibrium field place the magnetic axis at lower `Z` than the gfile axis, even when the LCFS shape itself remains non-symmetric?
