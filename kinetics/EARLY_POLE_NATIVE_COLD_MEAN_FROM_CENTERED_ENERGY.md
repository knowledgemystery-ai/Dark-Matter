# Stored centered energy resolves the sign of the native mean momentum

The [shared-bin analysis](EARLY_POLE_NATIVE_COLD_SHAPE_ATTRIBUTION.md) could resolve the occupation-only mean shift from panel masses and support, but that information alone allowed either sign for the full mean. The same complete numerical ledgers also retain a centered-energy moment per panel. It narrows the mean without reconstructing a profile inside the panel.

For the held scalar mass m write u=E−m=q²/(E+m). Given the panel's number mass and centered-energy mass, their ratio is its mean u. The function q(u)=sqrt(u(2m+u)) is concave. Its chord through the two native endpoints gives a lower bound on mean q; Jensen gives an upper bound. Both follow for any nonnegative intrapanel measure. Independently verify the original total-energy moment m times number plus centered energy. Remove the same finite Bose factor from both moments before normalizing the raw cases, and use the shared conditional panel shape when crossing the scalar histogram.

| D4 minus D2 mean q | Forward lower to upper | Inverse lower to upper |
|---|---:|---:|
| Raw | 2.70196e−7 to 2.72986e−7 | 2.64917e−7 to 2.67622e−7 |
| Native weighted | 3.64383e−8 to 5.62007e−8 | 4.03359e−8 to 6.00055e−8 |
| Explicit occupation-only contribution | 1.39622e−8 to 1.41024e−8 | 1.40072e−8 to 1.41444e−8 |
| Kernel-only contribution | 2.24015e−8 to 4.21729e−8 | 2.62561e−8 to 4.59337e−8 |

The positive native mean difference follows from the stored **moments**, and is not supported by panel counts alone. The two numerical meshes reproduce the interval endpoints within 1.2e−17 absolute. These are conditional held-chart means within the complete retained band, not a physical spectrum, a D8 result, or a spatial error estimate.

```bash
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error -m unittest test_early_pole_native_cold_mean_from_centered_energy
d=2026-11-28
p=results/early-pole-native-cold
python -W error early_pole_native_cold_mean_from_centered_energy.py results \
  "$p-shape-attribution-$d.json" | cmp - "$p-mean-from-centered-energy-$d.json"
```
