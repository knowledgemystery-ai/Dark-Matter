# Separating explicit Bose weights from the full cold-band rate contrast

The [two-point interventions](EARLY_POLE_NATIVE_COLD_CROSS_FACTORS.md) attribute a raw coefficient change to the shifted scalar pole at specified q. The [complete native cold ledgers](INDEPENDENT-EARLY-POLE-ENDPOINT-FINE-D4-COLD-CONTINUUM-LADDER-2026-11-28.md), by contrast, already multiply by the finite scalar occupation f. Use both accepted complete meshes on each chart. Every native q panel lies in one exact common scalar bin, so divide its forward integral by 1+f and its inverse integral by f. For the stimulated slope subtract the unweighted forward term before dividing by f. This exactly undoes the stored factors without modeling the intrapanel source shape, and retains the native q²/E phase measure. Then apply each archived scalar histogram to either kernel.

| Retained band | D4/D2 forward without explicit Bose factors | Native weighted D4/D2 | D4/D2 with the D2 histogram held on both | Explicit histogram-change factor, average over both orders |
|---|---:|---:|---:|---:|
| Full [1e−6,1.5e−5] | .09966014 | .05872260 | .06028798 | .97366567 |
| Lower [1e−6,1e−5] | .07737489 | .05872233 | .06028778 | .97366465 |
| Upper [1e−5,1.5e−5] | .10349978 | .11446238 | .10129792 | 1.12995657 |

The full raw integral is already about ten times smaller in D4. Applying even the **same** scalar histogram to both kernels makes the disparity larger, because it weights their different q profiles. Switching the histogram alone contributes another forward factor .9737 to the full comparison, but it **increases** the exact upper subband by about 13.0% in either order. Those statements separate the effect of a common weighting from the effect of changing the histogram. They are not an independent full-band pole attribution.

The finite affine slope (fD+A) is used throughout, not the large-f limit. The held diagnostic photon occupation from the ladder gives opposite full-band signs on D2 and D4 with **either** scalar histogram; both upper nets remain negative in both cases. The two meshes agree for all new unweighted and crossed sums within 1.73e−13. No actual photon occupation, new scalar trajectory, or spatial convergence follows.

```bash
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error -m unittest test_early_pole_native_cold_bose_reweight
python -W error early_pole_native_cold_bose_reweight.py results |
  cmp - results/early-pole-native-cold-bose-reweight-2026-11-28.json
```
