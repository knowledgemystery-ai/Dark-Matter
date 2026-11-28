# One held occupation cannot extend the full-band cold sign reversal to the upper band

The [earlier subband check](EARLY_POLE_ENDPOINT_FINE_D4_COLD_SUBBAND_ENCLOSURE.md) evaluated one common occupation. Instead calculate the neutral threshold for **every** constant held occupation from both accepted quadratures on D2 and D4. The conditional net is nD−B, with positive stimulated slope D and inverse B. For a subband with crossings inside native panels, the independent positive panel bounds give B_lower/D_upper ≤ n_neutral ≤ B_upper/D_lower. Take the union hull across the two observed meshes. The reversal D2 net negative and D4 net positive requires n strictly between the D4 and D2 thresholds. An interval between the inner endpoints certifies the reversal; an interval between the outer endpoints only says the bounds do not rule it out.

| Band q, in units of 1e−6 | Open held occupation interval certifying D2 loss and D4 gain, in units of 1e−9 |
|---|---:|
| Full retained [1,15] | (2.29937066, 2.90893020) |
| Central [3,5] | (2.45777300, 2.47366957) |
| Exact upper [10,15] | (2.92040842, 3.55121125) |

The full and central intervals overlap, although the original diagnostic occupation 2.5862538e−9 is not in that central guaranteed interval. The full and upper intervals do **not** overlap even using the outer, permissive bounds. The first upper neutral threshold sits 1.148e−11 above the latest possible end of the full-band reversal. Thus whenever the two accepted full-band summaries have opposite signs in this direction, **both** upper-band nets are negative for the same constant occupation. The other three tested subbands lack a certified interval of that sign reversal without additional intrapanel information.

This is a separation between two finite spatial charts and retained panel sums, not a spatial error bound, a constraint on an actual photon occupation, or a prediction for momentum-dependent live populations.

```bash
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error -m unittest test_early_pole_endpoint_cold_neutral_occupation_corridors
python -W error early_pole_endpoint_cold_neutral_occupation_corridors.py results |
  cmp - results/early-pole-cold-neutral-occupation-corridors-2026-11-28.json
```
