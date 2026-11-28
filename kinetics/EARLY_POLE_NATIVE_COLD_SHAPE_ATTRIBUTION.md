# The frozen kernel and explicit occupancy shifts reinforce on the shared cold bins

Use the [finite Bose cross](EARLY_POLE_NATIVE_COLD_BOSE_REWEIGHT.md) on both complete D2/D4 meshes. Normalize the two native spectra and the two crossed spectra independently on the complete held band. Compare their ten exact shared scalar bins. Average the signed kernel difference over both scalar histograms and the signed histogram difference over both kernels. These vectors close to the native D4-minus-D2 difference.

| Measure | Shared-bin total variation | Kernel contribution | Explicit occupancy contribution |
|---|---:|---:|---:|
| Forward number | 1.978563% | 1.157183 points | .821380 points |
| Inverse number | 2.086208% | 1.277152 points | .809057 points |
| Forward centered energy | 1.226421% | .728512 points | .497909 points |

The signed contributions agree bin by bin. The bin q=[2.28736e−6,3.30600e−6] accounts for 99.88% of lost normalized forward share. This concerns shape, not the roughly 17-fold decrease in total weighted rate. The kernel includes the held pole, charged state and momentum measure together.

For the occupation-only mean shift, the two histograms share each panel's conditional shape on a fixed kernel. Multiply the signed normalized weight difference by either exact panel endpoint to enclose its contribution. The forward mean shift is positive, in [1.08177e−8,1.71888e−8]; the full native mean shift may have either sign, in [−3.58222e−7,4.94417e−7]. These are support-only bounds on held numerical charts, not physical weighting, a D8 result, or trajectory attribution.

```bash
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error -m unittest test_early_pole_native_cold_shape_attribution
d=2026-11-28
python -W error early_pole_native_cold_shape_attribution.py results \
  "results/early-pole-native-cold-bose-reweight-$d.json" |
  cmp - "results/early-pole-native-cold-shape-attribution-$d.json"
```
