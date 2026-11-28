# The finite D4 daughter support clips only the upper inverse aperture

The [crossed raw point test](EARLY_POLE_NATIVE_COLD_CROSS_FACTORS.md) finds that the changed scalar pole dominates the D2-to-D4 reduction at two held momenta. Independently inspect the kinematic angular support imposed by the **finite retained numerical charts**. For forward absorption take the parent energy cap; for inverse take the daughter energy cap. Both last occupied cells are positive. The channels have different caps, about 32 and 26 in the archive's dimensionless units.

For the past-directed transfer t=E_scalar−ω, let E_p−E_d=t at the cap. The minimum momentum separation is p_d−p_p. It decreases with increasing allowed energy on this branch; the sum of momenta at the cap is much larger than the transferred spatial magnitude. Therefore the maximal allowed cosine is determined by the exact cap condition |k−q| crossing the separation at cosine +1. Across the retained q interval [1e−6,1.5e−5], the derivative of separation minus (k−q) stays above .9826 in all four tests. This proves whether clipping starts, not just whether it is visible at a sampled point.

| Native chart | Finite-support channel | First clipping momentum | Fraction of cosine interval retained at q=1.5e−5 |
|---|---|---:|---:|
| D2 | Forward parent | none in band | 1 |
| D2 | Inverse daughter | none in band | 1 |
| D4 | Forward parent | none in band | 1 |
| D4 | Inverse daughter | 1.2639838e−5 | .92329938 |

The D4 inverse loses about 7.67% of the **geometric cosine interval** at the top, not 7.67% of its rate. The cap rule is separately reproduced by the original energy-threshold solver at the limited angle. The two earlier point comparisons are below the clipping momentum, so their large raw suppression cannot be the loss of this part of the cosine interval. No conclusion follows about occupations beyond the numerical caps or their effect on a physical photon source.

```bash
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error -m unittest test_early_pole_native_cold_retained_aperture
d=2026-11-28
python -W error early_pole_native_cold_retained_aperture.py results \
  "results/early-pole-native-d2d4-cold-cross-q41-$d.json" |
  cmp - "results/early-pole-native-cold-retained-aperture-$d.json"
```
