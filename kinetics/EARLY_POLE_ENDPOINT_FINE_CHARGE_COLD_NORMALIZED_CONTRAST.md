# The cold-source magnitude collapse and the shape change are distinct

Compare the complete fine D1/D2/D4 native panel ledgers on the exact common q interval [1e−6,1.5e−5]. Track total response ratios separately from the three normalized positive measures: gross forward number, gross inverse number, and gross forward centered energy. For normalized profiles, allow arbitrary placement within each native panel but retain its known mass. No linear interpolation through panel interiors is used.

| Pair | Observable | Total-response ratio | Minimum compatible normalized TV | Minimum compatible W1 in q | Maximum compatible W1 in q |
|---|---|---:|---:|---:|---:|
| D1→D2 | Forward number | .201708 | .86508% | 1.311e−10 | 8.847e−7 |
| D1→D2 | Inverse number | .175191 | .91702% | 1.532e−10 | 8.879e−7 |
| D1→D2 | Forward centered energy | .204933 | .55083% | 4.712e−10 | 9.440e−7 |
| D2→D4 | Forward number | .058723 | 1.97856% | 2.167e−9 | 4.944e−7 |
| D2→D4 | Inverse number | .046417 | 2.08621% | 2.594e−9 | 4.979e−7 |
| D2→D4 | Forward centered energy | .060593 | 1.22642% | 5.807e−9 | 4.883e−7 |

TV is the minimum fraction that cannot be shared by any two measures satisfying their exact native interval capacities. W1 measures distance through q, rather than mass mismatch alone. Both its minimum and its maximum reflect permitted intrapanel placements rather than the true shape. Minima have independent primal/dual brackets; maxima have limiting endpoint-placement witnesses. Both separately accepted numerical quadratures give the same conclusions (TV differences ≤3.64e−14, transport differences ≤1.24e−19 in q).

Thus scalar normalization does **not** make the D2 and D4 shapes compatible, even though the main effect of charge refinement is the much larger change in total rate. Normalized TV bounds do not estimate error in total magnitude, and the wide W1 interval does not resolve the actual intrapanel displacement. Neither quantity is a spatial convergence or physical production claim.

```bash
cd kinetics
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
python -W error early_pole_endpoint_fine_charge_cold_normalized_contrast.py results \
  results/early-pole-fine-d4-cold-continuum-ladder-2026-11-28.json
```
