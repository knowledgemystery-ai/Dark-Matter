# The terminal amplitude-dependent sign does not require a particular cold cutoff

Use the exact five completed 4096-step trajectories from the [time-refinement study](EARLY_POLE_ENDPOINT_SYMMETRIC_CHARGE_TIME_REFINEMENT.md) and the [half-amplitude study](EARLY_POLE_ENDPOINT_SYMMETRIC_CHARGE_HALF_AMPLITUDE.md), without altering or rerunning any trajectory. Recover the scalar number at all 67 native edges from their terminal checkpoints, verifying the saved whole-spectrum count, the original cold count, and the pole in each run. Compare every perturbed cumulative count with the same nominal cumulative count.

| Cumulative terminal number versus nominal | Full lower | Full upper | Half lower | Half upper |
|---|---:|---:|---:|---:|
| q ≤ 1e−5 | −5.66913e−6 | −1.06050e−6 | −2.02814e−6 | +3.47856e−7 |
| Whole retained scalar chart | −5.66921e−6 | −1.06059e−6 | −2.02816e−6 | +3.47835e−7 |

All four signs match this pattern on the last 47 native edges, from q=2.28736e−6 up to the chart endpoint q=25. At the preceding edge q=1.58258e−6 the half-upper contrast is negative instead. These are exact retained cell-edge sums, not a statement about pointwise spectra.

The interior of a native cell need not have the same signs as its two edges. For an arbitrary cutoff inside a cell, take the difference at the lower edge, subtract at most that cell's nominal mass, and add at most that cell's perturbed mass. This uses only nonnegative allocations, with no common profile assumption. All four signs are then certified, conditional on the saved cell masses, for **every** threshold between 9.9817731491e−6 and 25. The least positive half-upper margin in these enclosures is 3.24597e−7. In the preceding open cell, 6.90620e−6 to 9.98177e−6, these broad bounds guarantee the two lower-branch signs only. Do not extend the certified range inward to the first positive edge. This is a cutoff robustness result for saved trajectories, not spatial convergence, a reconstruction of the distribution, or a physical production prediction.

```bash
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error -m unittest test_early_pole_symmetric_charge_cutoff_audit
python -W error early_pole_symmetric_charge_cutoff_audit.py results |
  cmp - <(gzip -dc results/early-pole-symmetric-charge-cutoff-audit-2026-11-28.json.gz)
```
