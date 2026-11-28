# A small momentum dependence breaks the constant-occupation obstruction

The [neutral threshold audit](EARLY_POLE_ENDPOINT_COLD_NEUTRAL_OCCUPATION_CORRIDORS.md) excludes any one constant held occupation from simultaneously making D2 lose and D4 gain on the full retained band and its exact upper subband. Both ledgers have an exact common native edge at q=1e−5, so split them there without fitting either panel. Use one common two-zone occupation in both grids: n_low below the edge and n_high above.

For **any** common bounded shape producing the stated signs on all four accepted ledgers, its maximum must exceed the largest D4 upper neutral threshold. Its minimum must be below the smallest D2 full-band neutral threshold. Positivity of the slopes therefore forces a strict peak-to-trough range greater than **1.14782256e−11**. This uses no assumption on how the shape varies within a panel.

Within the two-zone class the exact infimum is **1.14782816e−11**, only 4.88 ppm above that general lower bound on the required range. To see this, let n_high approach the largest D4 upper neutral threshold from above. At fixed n_high, D2 full-band loss imposes an upper bound on n_low, decreasing linearly with n_high in both meshes. The largest admissible n_low at the limiting n_high is therefore optimal; D4 full-band gain and D2 upper loss both have positive slack there. The limit is not attained because the inequalities are strict.

A constructive interior test uses n_low=2.908930039244118e−9 and n_high=2.920408520825516e−9. The four desired signs are reproduced in **both** meshes of each chart. The sharp lower-bound calculation is a diagnostic on fixed coefficient ledgers; a near-minimal construction has small sign margins. It is not an inferred shape, live radiation calculation, experimental constraint, or guarantee of spatial convergence.

```bash
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error -m unittest test_early_pole_endpoint_cold_two_zone_occupation
d=2026-11-28
p=results/early-pole-cold
python -W error early_pole_endpoint_cold_two_zone_occupation.py results \
  "$p-neutral-occupation-corridors-$d.json" | cmp - "$p-two-zone-occupation-$d.json"
```
