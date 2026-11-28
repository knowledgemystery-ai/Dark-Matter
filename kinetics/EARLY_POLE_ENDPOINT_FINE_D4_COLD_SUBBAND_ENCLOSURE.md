# D4 suppresses gross cold-band contributions in every tested subband

Use the complete accepted D2 and D4 fine meshes, not an interpolation over their different native panel boundaries. Enclose the weights in each subband by summing panels wholly inside for the lower bound and all intersecting panels for the upper bound. Positive gross forward, inverse, and slope ledgers allow a separate interval enclosure for each contribution and for the conditional net.

| q band (units of 1e−6) | D4/D2 gross forward enclosure | D4/D2 gross inverse enclosure | Conditional net D2 | Conditional net D4 |
|---:|---:|---:|:---:|:---:|
| [1, 3] | .03988–.38523 | .03121–.30390 | unresolved | unresolved |
| [3, 5] | .04957–.06237 | .03934–.04926 | unresolved | positive |
| [5, 8] | .06666–.08845 | .05347–.07045 | unresolved | unresolved |
| [8, 10] | .08926–.11531 | .07250–.09306 | negative | unresolved |
| [10, 15] | .11446238 | .09413039 | negative | negative |

The absolute forward and inverse reductions hold in all five tested subbands without assuming an intrapanel shape. The net signs use the same held photon occupation 2.5862538e−9 as the full-band D2/D4 midpoint test. Despite opposite full-band net signs, both charts are negative on the exact upper subband; the full-band positive D4 result is not a uniform-q sign result. “Unresolved” means only that the positive panel ledgers do not separate the interval from zero. No spatial convergence bound or physical occupation is inferred.

```bash
cd kinetics
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
python -W error -m unittest test_early_pole_endpoint_fine_d4_cold_subband_enclosure
python -W error early_pole_endpoint_fine_d4_cold_subband_enclosure.py \
  results/early-pole-fine-d4-cold-continuum-ladder-2026-11-28.json \
  results/fine-d2-endpoint-continuum-cold-2.jsonl.gz results/fine-d4-endpoint-continuum-cold-2.jsonl.gz
```
