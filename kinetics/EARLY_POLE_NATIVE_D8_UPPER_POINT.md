# A separate upper point checks the raw coefficient and real-Pauli arithmetic

Select q=1.3344923954e−5 before the calculation, using the midpoint of the shared D2/D4 raw median enclosure in the [raw-versus-weighted audit](EARLY_POLE_NATIVE_COLD_WEIGHT_TARGET.md). This is within the [planned upper band](EARLY_POLE_NATIVE_D8_RAW_BAND_PLAN.md), and disjoint from the active central D8 band. At k=.0015 on the exact independent x=.05 states, compute the complete angular average at this one q, with stock and real-Pauli-plus-mirror arithmetic at each of two angular meshes.

| Contrast | Raw gross forward | Raw gross inverse |
|---|---:|---:|
| D4/D2 | .10345413 | .08499090 |
| D8/D4 | .28292253 | .25536913 |

Relative agreement between implementations is within 2.89e−15 at the same mesh; the two observed meshes agree within 2.75e−13. The direct algorithm comparison agrees on angular support counts as well. The inverse geometric cosine fractions are 1, .97420 and .82599, in D2, D4 and D8 respectively. These are not the inverse rate ratios.

This is an actual direct raw angular evaluation at **one fixed q**, not a whole native panel or raw-band integral. No upper-band production job was started, no optional method was enabled for it, and the original central calculation was not altered. Do not use the point to infer band convergence, an actual photon population or a physical rate.

```bash
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error -m unittest test_early_pole_native_d8_upper_point
# This fresh comparison can take several minutes in a shared workspace.
python -W error early_pole_native_d8_upper_point.py results |
  cmp - results/early-pole-native-d8-upper-point-2026-11-28.json
```
