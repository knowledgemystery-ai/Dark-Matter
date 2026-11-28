# Replicate the real-Pauli arithmetic comparison at the first selected upper D8 panel

The [fine upper canary](EARLY_POLE_NATIVE_D8_UPPER_FINE_CANARY.md) used interior upper panel 183 at mesh ((8,10,16,24)). Use a **different** native panel, index 157, at mesh ((4,6,10,18)). This panel begins at the exact lower edge of the [separately planned upper raw span](EARLY_POLE_NATIVE_D8_RAW_BAND_PLAN.md), (q=1.0689130895889085	imes10^{-5}), and ends at (1.0787269717785436	imes10^{-5}). Bind the plan bytes, native D8 held state, full stock topology, pole, and real-Pauli implementation.

Compare original stock angular evaluation to the real-Pauli kernel with exact azimuth reflection. On this panel the largest difference over the nine integrated weighted coefficients is (4.00	imes10^{-15}), below the preregistered (2	imes10^{-12}) threshold. Replaying four archived real node values reproduces the complete receipt and rejects a checkpoint with changed plan metadata. Because the two canaries use different panels, this is **not** a same-panel angular-mesh convergence test or a sum over the planned 38 panels. The live central-band ledger is not read, changed, or replaced.

```bash
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error -m unittest test_early_pole_native_d8_upper_edge_coarse_check
d=2026-11-28; scratch=$(mktemp -d)
cp "results/d8-upper-edge-coarse-checkpoints-$d/"* "$scratch/"
python -W error early_pole_native_d8_upper_edge_coarse_check.py results "$scratch" |
  cmp - "results/early-pole-native-d8-upper-edge-coarse-check-$d.json"
```
