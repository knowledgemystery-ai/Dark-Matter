# The planned upper D8 inverse aperture is clipped throughout

The [retained D2/D4 map](EARLY_POLE_NATIVE_COLD_RETAINED_APERTURE.md) isolated the finite daughter cap in the inverse channel. Apply the same independent collinear test to the exact native D8 state used by the [disjoint raw plan](EARLY_POLE_NATIVE_D8_RAW_BAND_PLAN.md). No D8 collision integrals are evaluated here.

| Chart and diagnostic | D4 inverse | D8 inverse |
|---|---:|---:|
| First clipped q | 1.26398378e−5 | 8.60581691e−6 |
| Geometric cosine fraction at q=1.5e−5 | .92329938 | .79144336 |

The original running central D8 interval [3e−6,5e−6] lies entirely below the D8 clipping onset. The planned raw interval [1.06891309e−5,1.5e−5] lies entirely above it. The forward parent interval remains fully open throughout both. On the raw target the D8 inverse retains between **79.14% and 90.43%** of the geometric cosine interval (rounded). The end values bracket the whole interval: differentiating c=(k²+q²−v²)/(2kq), with v the cap momentum separation, gives a strictly negative upper bound on its derivative numerator. A separate bound keeps the inverse aperture positive.

These percentages are angular domain lengths, not reaction weights. They do not establish the D8 normalized raw profile or the impact of omitted occupation beyond the finite chart. No live central job was changed.

```bash
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error -m unittest test_early_pole_native_d8_retained_aperture
d=2026-11-28
python -W error early_pole_native_d8_retained_aperture.py results \
  "results/early-pole-native-d8-raw-band-plan-$d.json" |
  cmp - "results/early-pole-native-d8-retained-aperture-$d.json"
```
