# Locate the charge-grid cold-count difference in time and momentum

Use the 16 independently evolved native fine charge checkpoints (D1, D2, D4, D8 at x=.01, .02, .04, .05). The scalar grid and q≤1e−5 cold cut are identical. Compute the cold-count contrast between adjacent charge levels for every matched snapshot, then subtract successive contrasts. These signed increments telescope to the exact terminal contrast; they are not collision rates or time-discretization errors.

| Charge step | Already present at x=.01 | Added .01→.02 | Added .02→.04 | Added .04→.05 |
|---|---:|---:|---:|---:|
| D1→D2 | .00497% | .837% | 58.995% | 40.163% |
| D2→D4 | .00104% | .199% | 26.362% | 73.439% |
| D4→D8 | .00032% | .173% | 28.106% | 71.721% |

Each row is normalized to that pair's terminal signed cold-count deficit. The last two refinements each accumulate over 70% of that deficit in the final .04→.05 interval, unlike D1→D2. The two neighboring scalar cells spanning approximately q=[2.2874e−6,4.7783e−6] supply 102.75% and 105.78% of the terminal signed deficits for the latter two pairs; increases in other cells offset some of their decreases. This use of a signed denominator does not claim more than 100% of the absolute differences. For D4→D8 alone the lower of these two scalar cells accounts for 84.85% of the net terminal cold-count change.

Independent charged ledgers use each chart's own quadrature weights: at every interval the scalar and daughter number increments agree with the opposite parent increment to 6.60e−17 absolute or better. This localizes where the observed numerical discrepancy accumulated without assigning a cause to a particular collision term. It motivates a crossed terminal continuation, not an assertion that all earlier differences disappear if that interval alone is refined.

```bash
cd kinetics
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
python -W error early_pole_endpoint_charge_cold_count_history_budget.py results \
  results/early-pole-endpoint-charge-cold-count-gap-ladder-2026-11-28.json
```
