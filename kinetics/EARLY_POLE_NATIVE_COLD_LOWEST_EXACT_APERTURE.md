# Directly install the last unresolved low cut in native D2 and D4

The [updated adjacent band ledger](EARLY_POLE_NATIVE_COLD_PATCH_ADJACENT_CORRIDORS.md) still allows either sign in the lowest subband (q=[1,3]	imes10^{-6}). Reintegrate only the two children of the single panel containing (q=3	imes10^{-6}) for each native chart. No other panel, state, held occupation, or angular kernel is changed.

| Chart | Native parent ID | Original edges (units of (10^{-6})) | New child IDs |
|---:|---:|---:|---:|
| D2 | 7 | 2.4595998339, 3.1356591258 | 7, 8 |
| D4 | 18 | 2.9745740513, 3.1916537495 | 18, 19 |

For each mesh, require positivity, affine and centered-energy closure independently on each of the two children. Then sum them and recover every component of the original parent within (6	imes10^{-5}). Rejoin untouched original panels and apply the same componentwise tolerance to the paired complete retained subband. Report a held sign only from the observed per-panel two-mesh net box, and check that the resulting box fits the original looser allowed range. This is not a rigorous quadrature-error proof, a physical occupation, or a D8 result.

```bash
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
scratch=$(mktemp -d)
for mode in plan d2m1 d4m1 d2m2 d4m2; do
  python -W error early_pole_native_cold_lowest_exact_aperture.py results "$scratch" "$mode"
done
python -W error early_pole_native_cold_lowest_exact_aperture.py results "$scratch" finish
```

## Completed held comparison

At the same held (n=2.5862538037	imes10^{-9}), both meshes put D2 below zero (net (-2.10557	imes10^{-17})) and D4 above zero (net (+1.52408	imes10^{-18})). The D2 and D4 band neutral values are (2.84587	imes10^{-9}) and (2.23070	imes10^{-9}), respectively. All newly inserted components close their original parents within (1.45	imes10^{-12}), and retained components differ across meshes by at most (4.19	imes10^{-13}). Both signed boxes remain inside their earlier loose numerical bounds.

The new signs fit the previous common test-occupation corridor ((2.43556,2.60358)	imes10^{-9}) without narrowing it. This is the same constant test occupation, not a recovered physical population.

```bash
d=2026-11-28
python -W error -m unittest test_early_pole_native_cold_lowest_exact_aperture
python -W error early_pole_native_cold_lowest_exact_aperture.py results \
  "results/cold-lowest-exact-aperture-patches-$d" finish |
  cmp - "results/early-pole-native-cold-lowest-exact-aperture-$d.json"
```
