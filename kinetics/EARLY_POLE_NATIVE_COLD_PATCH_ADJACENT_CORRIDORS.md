# The adjacent D4 [8,10] band loses at the original held occupation

The [D4 exact-split calculation](EARLY_POLE_NATIVE_D4_CLIPPED_NEXT_BAND.md) installed cuts at (5) and (8), in units of (10^{-6}). The original native mesh already has an edge at (10). Rejoin the exact archived D4 patch on both meshes and use the resulting exact ([8,10]) span as well. Do not use interpolated numbers from the old intersecting panels. D2 comparisons retain the original centered-moment enclosures.

| q band (units of (10^{-6})) | New D4 neutral occupancy (units of (10^{-9})) | Held D2 / D4 sign |
|---:|---:|:---:|
| [1,3] | 1.781–2.797 | unresolved / unresolved |
| [3,5] | 2.22969–2.38344 | negative / positive |
| [5,8] | 2.43556 | negative / positive |
| [8,10] | 2.64867 | negative / **negative** |
| [10,15] | 2.92041 | negative / negative |

At (n=2.5862538037	imes10^{-9}), the adjacent D4 ([8,10]) net is (-5.17230	imes10^{-22}) on the observed mesh box. The earlier centered-moment enclosure alone could not settle its sign. D4's **band-averaged** neutral occupation is ordered and separated across the four bands starting at 3; this does not establish a monotone pointwise threshold. Enclosing inverse and slope independently across both meshes, the entire set of resolved D2/D4 signs, including the full-band signs, holds over the sufficient open interval (nin(2.43556254,2.60358266)	imes10^{-9}). The lower limit comes from D4 ([5,8]); the upper from the D2 centered bound on ([3,5]). Both charts remain unresolved in the lowest subband with the saved data. The corridor is a held numerical diagnostic, not a live occupation, D8 result, or quantified continuum error.

```bash
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error -m unittest test_early_pole_native_cold_patch_adjacent_corridors
python -W error early_pole_native_cold_patch_adjacent_corridors.py results |
  cmp - results/early-pole-native-cold-patch-adjacent-corridors-2026-11-28.json
```
