# What changes the sign in the three controlled G9 pressure probes?

Replay one versus two local substeps at the three exact held starts. Use each previously recorded arm's actual endpoint occupation and the same common x=.05 mass at both measurements. Compare the published warm count-and-energy pressure reconstruction with a diagnostic based solely on native nodal weights: the sum of n_i q_i² / hypot(q_i,m). The nodal diagnostic applies at every node, including those representing unsplit cold cells, and is not a replacement for the reconstructed observable.

| Pulse inserted at x | After pulse: reconstructed / nodal | At x=.05: reconstructed / nodal | Endpoint count component | Endpoint within-cell component |
|:---|---:|---:|---:|---:|
| .04 | +6.3263e-10 / +6.3006e-10 | −2.4517e-10 / −2.4520e-10 | −2.4895e-10 | +3.7763e-12 |
| .04598633 | +6.0340e-10 / +6.0096e-10 | −2.7011e-10 / −2.7031e-10 | −2.6492e-10 | −5.1877e-12 |
| .04862462 | +1.3198e-9 / +1.3143e-9 | −1.0758e-10 / −1.0788e-10 | −1.0915e-10 | +1.5737e-12 |

After the pulse and at the endpoint the two pressure diagnostics agree to within 0.5% of the reconstructed contrast in all six comparisons. The cell-count and within-cell columns use the exact symmetric product rule on reconstructed cell pressure: each cell's pressure change is the average pressure per particle times the change in count plus the average count times the change in pressure per particle. Terms are normalized to that arm pair's one-step total pressure. At the endpoint the magnitude of the within-cell term is **at most 2.0%** of the count term. Cold cells still make a positive endpoint contribution; the negative pressure term from changed warm-cell populations outweighs it. This is an attribution of these discrete contrasts, not a causal identification of a specific reaction, an extrapolation to other clocks, or a physical convergence statement.

## Reproduction

From the kinetics directory, with the three existing published pulse receipts and corresponding checkpoints in results:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 python -W error \
  early_pole_endpoint_warm_propagated_closure_contrast.py results
```
