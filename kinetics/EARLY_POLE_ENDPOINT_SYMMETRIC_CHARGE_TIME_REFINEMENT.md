# Halving the time step preserves the sampled odd/even charge-shape behavior

Reintegrate the nominal charged remap and both truly equal-opposite moment-null shapes from the exact original .04 states. All three runs use the same terminal clock with 4096 steps instead of 2048. The scalar seed is unchanged. Both grids provide eight exact common saved times ending at .05, with attached full final states.

| Observable | 2048 steps | 4096 steps |
|---|---:|---:|
| Terminal cold half upper-minus-lower count | +2.3048923e−6 | +2.3043126e−6 |
| Terminal cold mean-minus-nominal count | −3.3650830e−6 | −3.3648174e−6 |

All 16 sampled odd/even total-yield signs agree and all three sampled sign-reversal brackets of the odd contrast are unchanged. At the common clocks, the maximum absolute inter-grid changes are 5.35e−11 for the odd yield response and 3.00e−11 for the even response. The terminal cold odd magnitude changes by about .0252%; the even magnitude changes by .0079%. Both terminal shape-control counts remain below their own same-clock nominal state, and the even yield response stays negative and decreases at each of the eight sampled times.

This is an actual same-start time-refinement check, not a continuum time-error bound, a count of unsampled zero crossings, a physical rate result, or a spatial-convergence statement.

```bash
cd kinetics
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
p=results/early-pole
d=2026-11-28
for branch in lower nominal upper; do
  gzip -dc "$p-charge-symmetric-time-$branch-x05-checkpoint-$d.json.gz" > "/tmp/symmetric-time-$branch.json"
done
python -W error early_pole_endpoint_symmetric_charge_time_comparison.py results \
 "$p-charge-symmetric-null-terminal-$d.json" "$p-symmetric-charge-history-budget-$d.json" \
 --records "$p-charge-symmetric-time-lower-x05-receipt-$d.json" "$p-charge-symmetric-time-nominal-x05-receipt-$d.json" "$p-charge-symmetric-time-upper-x05-receipt-$d.json" \
 --checkpoints /tmp/symmetric-time-{lower,nominal,upper}.json
```

To recompute each fine continuation, use `early_pole_endpoint_symmetric_charge_time_refinement.py results SUMMARY SIGN CHECKPOINT` with SIGN −1, 0, or 1 and the exact archived symmetric terminal summary.
