# Opposite moment-null inputs give a larger even than odd terminal cold response

The previous control moved the charged middle-node allocation by equal *fractions of two different feasible distances*. Here use a single signed direction in each common D2 coarse cell and **the same absolute displacement in both directions**: 5% of the smaller of the two distances from the nominal D2→D4 remap to its Pauli-feasible limits. Match both charge number and centered energy cell by cell, keep both shared middle occupations equal, and retain the exact D2 scalar seed. Continue the two affected states with the same 2048-step terminal rule.

| Observable | Lower | Nominal | Upper | Half upper-minus-lower | Mean endpoints minus nominal |
|---|---:|---:|---:|---:|---:|
| Frozen initial net scalar log-rate | 4.20482439 | 4.11264443 | 4.01507686 | −.09487377 | −.00269381 |
| Terminal cold count | .01705708349 | .01706275346 | .01706169327 | +2.30489e−6 | −3.36508e−6 |

Both final counts are again below nominal, even for a truly symmetric finite perturbation. The odd contrast switches sign from the frozen initial rate to the terminal count, and the terminal even contrast is larger in magnitude than the terminal odd one. That contrast concerns this finite amplitude, not a derivative, a rate decomposition, or all possible shapes. The larger endpoint difference is 1.7424% of the nominal charge-grid-plus-representation contrast. This supports sensitivity testing of the grid result, not a physical convergence claim.

```bash
cd kinetics
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
d=2026-11-28
p=results/early-pole
python -W error -m unittest test_early_pole_endpoint_charge_symmetric_null_terminal
for s in lower upper; do
  gzip -dc "$p-charge-symmetric-null-$s-x05-checkpoint-$d.json.gz" >"/tmp/symmetric-null-$s.json"
done
python -W error early_pole_endpoint_charge_symmetric_null_receipt.py results \
  "$p-endpoint-charge-grid-donor-terminal-$d.json" \
  --records "$p-charge-symmetric-null-lower-x05-receipt-$d.json" "$p-charge-symmetric-null-upper-x05-receipt-$d.json" \
  --payloads /tmp/symmetric-null-{lower,upper}.json
```

Use `early_pole_endpoint_charge_symmetric_null_terminal.py results SIGN CHECKPOINT`, with SIGN −1 or 1, to repeat the corresponding continuation.
