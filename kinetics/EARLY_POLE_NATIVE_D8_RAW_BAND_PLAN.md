# A separate upper-band plan for a future D8 raw-kernel check

The [central D8 run](EARLY_POLE_ENDPOINT_FINE_D8_CENTRAL_COLD_BAND_PROTOCOL.md) continues unchanged. The [raw-versus-weighted audit](EARLY_POLE_NATIVE_COLD_WEIGHT_TARGET.md) shows that its central band selects very little of the raw forward measure on the prior D2 and D4 charts. If raw coverage is needed later, use a disjoint upper target.

Using only the exact D8 full-band topology, enumerate every contiguous interval whose endpoints are native D8 edges. For each interval require at least 80% of **both** raw forward and raw inverse on **both** completed D2/D4 angular meshes. Count only whole contained source panels, so no unobserved within-panel profile is assumed. First minimize the number of D8 panels, then interval width and starting index. All selected D8 panels also lie in single native scalar bins, which permits a later exact Bose unweighting.

| Plan property | Value |
|---|---:|
| Candidate span q | [1.06891309e−5,1.5e−5] |
| Native zero-based D8 panels | 157 through 194, 38 of 195 |
| Candidate intervals searched | 19,110 |
| Guaranteed D2 raw forward / inverse coverage | at least 80.09% / 81.76% |
| Guaranteed D4 raw forward / inverse coverage | at least 84.11% / 85.67% |

The optimum is for the held topology and conservative lower sums only. These percentages are not D8 coverage estimates. **No new D8 rate was calculated and no active central job was altered.** Any later run would need its own exact paired-mesh complete check on the selected indices, with actual D8 raw totals reported only for that retained span.

```bash
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error -m unittest test_early_pole_native_d8_raw_band_plan
d=2026-11-28
python -W error early_pole_native_d8_raw_band_plan.py results \
  "results/early-pole-native-cold-weight-target-$d.json" |
  cmp - "results/early-pole-native-d8-raw-band-plan-$d.json"
```

## Future result admission

When warranted, the stock full-span pilot can evaluate just these planned IDs through its existing `panel_indices` argument in two separate output files. The original central-band ledger is not an input. The selected-band reducer binds the plan and reconstructed full topology, rejects incomplete, duplicated, mixed-mesh, or central-band rows, independently checks the affine and energy ledgers, and removes the finite scalar weight only after all expected rows are present. Its selected raw totals remain restricted to that band.

```bash
python -W error -m unittest test_early_pole_native_d8_raw_band_receipt
# Run this only after both future selected stock quadratures are complete.
python -W error early_pole_native_d8_raw_band_receipt.py results \
  "results/early-pole-native-d8-raw-band-plan-$d.json" \
  "$RAW_UPPER_STOCK_MESH_1" "$RAW_UPPER_STOCK_MESH_2"
```
