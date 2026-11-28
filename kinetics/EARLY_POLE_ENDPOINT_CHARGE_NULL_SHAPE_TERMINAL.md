# Two moment-preserving charge-shape alternatives leave the terminal grid contrast dominant

On the exact D2-to-D4 charge remap with D2 scalar seed at x=.04, move the charged allocation at the duplicated middle-energy nodes by ±5% of the feasible displacement from the nominal fit toward its corresponding endpoint. Redistribute outer-node masses analytically to preserve count and centered energy in every shared coarse cell. The two middle occupations remain equal and all Pauli caps hold. The largest cellwise moment residual is 2.22e−16. The initial scalar root and total energy are unchanged.

| Charge-shape move | Initial net scalar log-rate | Terminal cold count | Change from nominal terminal remap |
|---:|---:|---:|---:|
| −5% | 4.21522199 | .01705678638 | −5.96709e−6 |
| Nominal | 4.11264443 | .01706275346 | 0 |
| +5% | 4.01507686 | .01706169327 | −1.06019e−6 |

The initial frozen rate moves in opposite directions, but after two independent complete 2048-step continuations both terminal counts are lower than nominal. The larger absolute terminal change is 1.834% of the nominal path-averaged charge-grid replacement effect, so this specific pair of shape changes does not overturn the grid-dominance finding. It is **not** a bound for arbitrary feasible remaps or a physical spectrum. The independently checked final occupations, continuity ledgers, static two-order comparison and verified remap are retained.

```bash
cd kinetics
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
p=results/early-pole
d=2026-11-28
python -W error -m unittest test_early_pole_endpoint_charge_null_face_clock
python -W error early_pole_endpoint_charge_null_face_clock.py results "$p-two-moment-charge-from2-x04-setup-$d.json"
for branch in lower upper; do
  gzip -dc "$p-charge-null-$branch-x05-checkpoint-$d.json.gz" > "/tmp/null-shape-$branch.json"
done
python -W error early_pole_endpoint_null_shape_terminal_receipt.py results \
  "$p-charge-null-face-clock-$d.json" "$p-endpoint-charge-grid-donor-terminal-$d.json" \
  --records "$p-charge-null-lower-x05-receipt-$d.json" "$p-charge-null-upper-x05-receipt-$d.json" \
  --payloads /tmp/null-shape-{lower,upper}.json
```
