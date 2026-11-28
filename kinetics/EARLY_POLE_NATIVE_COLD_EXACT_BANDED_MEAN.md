# Most of the D4 mean-q increase comes from shifting weight between exact bands

Use both completed [matching D2/D4 partitions](EARLY_POLE_NATIVE_COLD_EXACT_PAIRED_BANDS.md). In every panel the positive mass and stored centered-energy moment bound the mean of (q=sqrt{u(u+2m)}) between its chord and Jensen bounds. Use the exact scalar bin to remove the finite Bose factors for the raw version, and preserve forward and inverse as separate positive measures. Do not apply this construction to the signed net.

For a given band (j), let (p_j) be its probability and (mu_j) its mean q. The whole D4–D2 difference is the sum of two terms: (sum_j(p_{4j}-p_{2j})(mu_{4j}+mu_{2j})/2) for shifting weight between the five bands, and (sum_j(p_{4j}+p_{2j})(mu_{4j}-mu_{2j})/2) for changing conditional means within them. This symmetric identity compares the two held results; it does not identify causes.

| Positive measure | Between-band term, (10^{-9}) | Within-band term, (10^{-9}) |
|---|---:|---:|
| Weighted forward | 28.1257–28.3062 | 8.70582–27.0132 |
| Weighted inverse | 30.3984–30.5856 | 10.3193–28.5587 |
| Raw forward | 187.008–187.151 | 83.1529–85.8601 |
| Raw inverse | 177.750–177.884 | 87.1326–89.7638 |

The **between-band term exceeds the within-band term** even with their conservative intervals and both meshes. D4's conditional raw mean exceeds D2's in all five bands for both gross measures. For the weighted means the four bands above (3	imes10^{-6}) also resolve upward, but the lowest band does not: the forward difference lies between (-2.92) and (9.85) in units of (10^{-9}). The whole-band upward result is therefore not an assertion that every weighted band mean increases. The allocation depends on the chosen band boundaries, not on a new dynamical calculation.

```bash
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error -m unittest test_early_pole_native_cold_exact_banded_mean
d=2026-11-28
python -W error early_pole_native_cold_exact_banded_mean.py results |
  cmp - "results/early-pole-native-cold-exact-banded-mean-$d.json"
```
