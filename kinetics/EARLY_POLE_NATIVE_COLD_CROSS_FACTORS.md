# The shifted scalar pole dominates two D2-to-D4 raw cold-rate differences

The [complete frozen cold integrals](INDEPENDENT-EARLY-POLE-ENDPOINT-FINE-D4-COLD-CONTINUUM-LADDER-2026-11-28.md) suppress the D4 gross forward source relative to D2. To isolate mechanisms at two specified points, take exact independently evolved D2 and D4 states at x=.05. Recompute the complete transverse angular average for all eight choices of scalar mass, photon energy and residue, and reconstructed charged distributions. Use the same native cell occupations and exact-mirror real-Pauli integrator. Verify both angular settings independently for every row.

| q | Native D4/D2 forward | Scalar-pole factor | Photon-pole factor | Charged-distribution factor |
|---:|---:|---:|---:|---:|
| 4.1e−6 | .06019539 | .06024615 | .99915753 | .999999998 |
| 8.2e−6 | .07726356 | .07732313 | .99922952 | .999999997 |

Factors exponentiate the mean contribution to log coefficient over all replacement orders; their product reproduces the exact native ratio. The scalar-pole factor accounts for about 99.97% of the forward log suppression at both points. Inverse is likewise dominated by the scalar pole, with total ratios .04757420 and .06243808; the direct charged factors for inverse are 1.00019661 and .99990573. The parent threshold at cosine zero increases from about 18.12 to 20.60 at the central point. The D4 all-native coarse result at q=4.1e−6 exactly repeats the independently archived real-Pauli fixture. The largest observed two-mesh relative disagreement over all 16 pairs is below 4e−13.

This changes the **raw** coefficients at fixed q. It does not contradict a strong effect of the charged **history** on the scalar pole, provide a factorization of the entire cold band, determine a live Bose-weighted source, or establish spatial convergence. Crossed inputs are mathematical interventions, not separately self-consistent time trajectories.

```bash
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error -m unittest test_early_pole_native_cold_cross_factors
d=2026-11-28
for point in 41 82; do
  if [ "$point" = 41 ]; then q=4.1e-6; else q=8.2e-6; fi
  scratch="$(mktemp -d)"
  python -W error early_pole_native_cold_cross_factors.py results "$q" "$scratch" |
    cmp - "results/early-pole-native-d2d4-cold-cross-q$point-$d.json"
done
```
