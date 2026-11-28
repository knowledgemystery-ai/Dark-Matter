# Centered-energy constraints tighten the cold momentum transport distance

For each complete D2/D4 forward or inverse measure, use the stored centered-energy moment to bound CDF values at any interior cut as in the [partial-panel check](EARLY_POLE_NATIVE_COLD_SUBBAND_FROM_CENTERED_ENERGY.md). On an open interval between two cuts, monotonicity places either CDF between its left lower and right upper enclosure, including cuts that coincide with native edges by one-sided panel ownership. Integrate the resulting pointwise distance intervals with their exact widths. This encloses the one-dimensional Wasserstein-1 distance without fitting a density or assuming straight CDF segments. The [independent mean difference](EARLY_POLE_NATIVE_COLD_MEAN_FROM_CENTERED_ENERGY.md) supplies a stronger lower bound where applicable.

| Normalized held measure | Moment-constrained D2/D4 W1 bracket (q units) |
|---|---:|
| Weighted forward | [3.64383e−8,3.80950e−7] |
| Weighted inverse | [4.03359e−8,3.82058e−7] |
| Raw forward | [2.70196e−7,4.99545e−7] |
| Raw inverse | [2.64917e−7,4.93115e−7] |

For the weighted forward comparison, the support-only bracket was about [2.167e−9,4.944e−7]: the new lower bound is about 16.8 times larger, and the upper bound about 23% smaller. The raw lower bounds are much larger, but the raw and weighted brackets still overlap; they do not order the true distances. Nested cut grids tighten the CDF contribution and both accepted angular meshes agree. These are conservative brackets, not an exact moment-constrained optimizer, a physical trajectory or a convergence claim.

```bash
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error -m unittest test_early_pole_native_cold_centered_energy_transport
d=2026-11-28
p=results/early-pole-native-cold
python -W error early_pole_native_cold_centered_energy_transport.py results \
 "$p-mean-from-centered-energy-$d.json" | cmp - "$p-centered-energy-transport-$d.json"
```
