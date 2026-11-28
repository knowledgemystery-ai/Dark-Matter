# Match the D2 band geometry to the exact five-band D4 partition

The [previous completed D4 geometry](EARLY_POLE_NATIVE_D4_EXACT_PARTITION.md) has all five D4 bands exact but the middle D2 bounds still allow intrapanel redistribution. Keep the already evaluated D2 cut at (3	imes10^{-6}). Install two further edges at (5	imes10^{-6}) and (8	imes10^{-6}), which intersect original D2 parent IDs 16 and 25. The four recalculated children of the final 52-panel native geometry are IDs 17, 18, 27, 28. The native edge at (10	imes10^{-6}) remains unchanged. Reuse all other archived components.

Require both children to reproduce each parent's stored number, energy and centered-energy columns separately, including each small inverse term. Require the observed componentwise tolerance (6	imes10^{-5}) for each retained paired mesh entry, not just whole-band totals. Exact partition closure must count every panel once. The common diagnostic occupation is held and is not a live photon population. Two observed quadratures do not bound global numerical or spatial error.

```bash
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
scratch=$(mktemp -d)
python -W error early_pole_native_d2_exact_partition.py results "$scratch" plan
python -W error early_pole_native_d2_exact_partition.py results "$scratch" mesh1
python -W error early_pole_native_d2_exact_partition.py results "$scratch" mesh2
python -W error early_pole_native_d2_exact_partition.py results "$scratch" finish
```

## Completed paired D2 quadrature

The native counts are 8, 10, 10, 7 and 17. All five D2 band nets are negative at the original held diagnostic. The central ([3,5]	imes10^{-6}) band supplies 86.35% of the full D2 loss magnitude. Both recomputed parents close within (1.14	imes10^{-13}); the largest observed per-panel component difference between the two meshes is (1.37	imes10^{-11}). See the [matched D2/D4 comparison](EARLY_POLE_NATIVE_COLD_EXACT_PAIRED_BANDS.md) for the exact neutral thresholds.

```bash
d=2026-11-28
python -W error -m unittest test_early_pole_native_d2_exact_partition
python -W error early_pole_native_d2_exact_partition.py results \
  "results/d2-exact-partition-patches-$d" finish |
  cmp - "results/early-pole-native-d2-exact-partition-$d.json"
```
