# Halving the charge-shape amplitude changes which cold endpoint lies above nominal

Keep the exact D2 scalar seed, the tested moment-null charged direction on the D4 destination grid, and the original 4096-step terminal clock. Rerun the two equal-opposite signs at half their original amplitude. All state differences, time samples, and scalar outcomes are compared to the same independently archived 4096 nominal trajectory.

| Terminal cold contrast | Full tested amplitude | Half amplitude | Ratio to linear/quad scaling from full |
|---|---:|---:|---:|
| Odd, half upper-minus-lower | +2.30431257e−6 | +1.18799749e−6 | 1.03111 of half |
| Even, mean minus nominal | −3.36481740e−6 | −8.40141295e−7 | .998736 of quarter |

At full amplitude both terminal cold counts fell below nominal; at half amplitude the lower arm is below by 2.02814e−6 while the upper arm is **above** by 3.47856e−7. The negative even term still exists, but shrinks enough that the odd term is larger. On the eight saved total-yield clocks, the even term's half-amplitude comparison lies within 0.162% of quadratic scaling; the odd comparison has up to 5.22% departure from linear scaling. All three sampled odd sign-reversal brackets are retained. These are two finite amplitudes on one numerical history, not infinitesimal derivatives, an asymptotic Taylor bound, spatial convergence, or a physical prediction. The direction and numerical time mesh have not changed.

```bash
cd kinetics
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
d=2026-11-28
p=results/early-pole
python -W error -m unittest test_early_pole_endpoint_symmetric_charge_amplitude_receipt
for branch in lower upper; do
  gzip -dc "$p-symmetric-half-$branch-x05-checkpoint-$d.json.gz" > "/tmp/symmetric-half-$branch.json"
done
python -W error early_pole_endpoint_symmetric_charge_amplitude_receipt.py results \
 "$p-charge-symmetric-time-comparison-$d.json" "$p-charge-symmetric-time-nominal-x05-receipt-$d.json" \
 --records "$p-symmetric-half-lower-x05-receipt-$d.json" "$p-symmetric-half-upper-x05-receipt-$d.json" \
 --checkpoints /tmp/symmetric-half-{lower,upper}.json
```

The corresponding two 4096-step trajectories can be replayed with `early_pole_endpoint_symmetric_charge_half_amplitude.py results SUMMARY SIGN CHECKPOINT`.
