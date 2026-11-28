# The central cold band targets the Bose-weighted rate, not the raw kernel

Use the [same finite Bose reweighting](EARLY_POLE_NATIVE_COLD_BOSE_REWEIGHT.md) on both complete D2/D4 cold quadratures. Normalize the positive forward number measure separately before and after the factor 1+f. Inner and outer sums of whole panel weights bound any partial subband. No source shape is assumed inside a panel.

| Chart | Raw forward share in q=[3e−6,5e−6] | Weighted share in that band | Raw share in exact upper [1e−5,1.5e−5] | Weighted share in exact upper |
|---|---:|---:|---:|---:|
| D2 | 0.726%–1.276% | 81.09%–95.65% | 85.30% | 4.88 ppm |
| D4 | 0.406%–0.572% | 80.75%–86.14% | 88.59% | 9.50 ppm |

The raw median lies in [1.304,1.346]e−5 for D2 and [1.325,1.344]e−5 for D4; the weighted median lies in [3.305,4.303]e−6 and [3.959,4.242]e−6. Since the Bose multiplier is constant on each native panel, the normalized raw-versus-weighted total variation within a chart is exactly the half-sum of the absolute normalized panel differences: .96884 in D2 and .97790 in D4. Between D2 and D4, the minimum compatible normalized shape difference is 4.026% for the raw measures versus 1.979% for the weighted measures, with primal and dual overlap certificates and two matching meshes.

The running [D8 central protocol](EARLY_POLE_ENDPOINT_FINE_D8_CENTRAL_COLD_BAND_PROTOCOL.md) remains a useful targeted **weighted** result. These two prior charts do not support interpreting it as a check of most of the raw kernel, making a D8 raw shape prediction, or replacing the original gate. Any later raw-focused quadrature would need its own held span and decision criteria.

```bash
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error -m unittest test_early_pole_native_cold_weight_target
d=2026-11-28
python -W error early_pole_native_cold_weight_target.py results \
  "results/early-pole-native-cold-bose-reweight-$d.json" |
  cmp - "results/early-pole-native-cold-weight-target-$d.json"
```
