# Check the G9 nodal pressure sign over a continuous diagnostic mass range

Use the six previously reproduced pulse and terminal nodal differences. Keep the saved occupations fixed and let only the scalar mass in the pressure weight vary continuously between the independently archived x=.04 and x=.05 native measurement masses, 0.0006770053032736531 and 0.0008731829968237441. This varies the observable, not the evolution or the pressure reconstruction.

Each nodal weight q²/hypot(q,m) decreases with m. Bound positive and negative contributions separately at the two mass endpoints, and sum with an explicit floating-point pad relative to the sum of absolute terms. The resulting single interval separates the sign from zero in each case, without sampling interior masses. Normalize all bounds for each comparison by the positive one-step nodal pressure at the original x=.05 reference mass.

| Insertion x | Bound immediately after the pulse | Bound at x=.05 |
|---|---:|---:|
| .04000000 | [+6.3005e−10, +7.6196e−10] | [−2.4573e−10, −2.4333e−10] |
| .04598633 | [+6.0095e−10, +7.2644e−10] | [−2.7032e−10, −2.6004e−10] |
| .04862462 | [+1.3143e−9, +1.5934e−9] | [−1.0912e−10, −9.1499e−11] |

These conditional nodal signs cannot be attributed to the particular common mass weight inside this interval. The enclosure does not bound evolution errors, a recalculated trajectory with a different mass, or the reconstructed pressure at the alternate masses.

```bash
cd kinetics
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
python -W error -m unittest test_early_pole_endpoint_warm_pulse_mass_sign_bound
python -W error early_pole_endpoint_warm_pulse_mass_sign_bound.py results
```
