# Compare the D2 and D4 cold nets on five identical exact bands

The [D2 exact partition](EARLY_POLE_NATIVE_D2_EXACT_PARTITION.md) and [D4 exact partition](EARLY_POLE_NATIVE_D4_EXACT_PARTITION.md) now share all six edges. Both sets of quadratures re-integrate the necessary native cut panels. No remaining mixture or distribution assumption about a panel intersecting one of these five bands is needed to compute its numerical net.

| q (units of (10^{-6})) | Neutral D2 | Neutral D4 | Signs at original held (n), D2/D4 |
|---:|---:|---:|:---:|
| [1,3] | 2.845873 | 2.230698 | −/+ |
| [3,5] | 2.915364 | 2.304691 | −/+ |
| [5,8] | 3.041886 | 2.435563 | −/+ |
| [8,10] | 3.262898 | 2.648675 | −/− |
| [10,15] | 3.551211 | 2.920408 | −/− |

Neutral occupations are shown in units of (10^{-9}). D4's neutral is below D2's in every matched band, although both lose in the upper two at the original held occupation. The observed full sign pattern holds throughout the open constant-occupation interval ((2.43556254,2.64867453)	imes10^{-9}). The stricter upper bound from the earlier loose D2 central interval no longer applies; both new limiting bands belong to D4, immediately adjacent at (q=8	imes10^{-6}). There is room for about 5.83% downward or 2.41% upward from the original diagnostic before one of these band signs becomes unresolved. These are observed two-mesh boxes on the exact held numerical bands, not continuum-error control, spatial convergence, or a determination of live photon occupation.

```bash
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error -m unittest test_early_pole_native_cold_exact_paired_bands
d=2026-11-28
python -W error early_pole_native_cold_exact_paired_bands.py results |
  cmp - "results/early-pole-native-cold-exact-paired-bands-$d.json"
```
