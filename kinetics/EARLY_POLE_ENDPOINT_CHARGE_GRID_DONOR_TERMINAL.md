# The destination charge grid dominates this fixed-scalar restart comparison

Repeat the .04→.05 terminal evolution with the D2 scalar occupation held **only at the start** in all four arms. Two arcs already exist: native D2 and native D4 charged grid/occupations with that scalar seed. For the remaining arcs, remap D2 charges to the D4 grid, and D4 charges to the D2 grid. Use one common D2 cell at a time. Match charge count and its centered energy moment at the actual restart mass, with no Pauli violations. The coarse two-endpoint remap is unique at these moments; the fine four-node remap uses an interior logit fit and is not a unique physical reconstruction. The total scalar root is unchanged and the largest relative cell moment discrepancy is 3.97e−14.

| Destination charge grid | Inherited D2 charge state | Inherited D4 charge state |
|---|---:|---:|
| D2 | .01738683682 | .01730606290 |
| D4 | .01706275346 | .01697933293 |

At x=.05 the cold-count difference on the diagonal of this **fixed-initial-scalar** square is −.0004075038862. The path-averaged grid-plus-representation change is −.0003254066626 (79.85%); the inherited charged-state change is −.0000820972237 (20.15%). The finite interaction is −.0000026466182 and is recorded separately from the path averages. Both effects are negative on both paths of the square.

This apportions the charge-only contrast conditional on one scalar restart and this chosen remap. It is not 79.85% of the full native D2→D4 difference, does not uniquely separate the grid from its discrete occupation representation, and does not establish robustness for arbitrary higher-moment remaps or physical convergence. The complete final occupations and conservative ledgers are retained.

```bash
cd kinetics
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
p=results/early-pole
d=2026-11-28
python -W error -m unittest test_early_pole_endpoint_charge_nested_moment_remap
for c in 2 4; do
  python -W error early_pole_endpoint_charge_nested_moment_remap.py results "$c" "/tmp/remap-checkpoint-from$c.json" --setup-only
  python -W error early_pole_endpoint_charge_nested_moment_remap.py results "$c" "/tmp/remap-checkpoint-from$c.json"
  gzip -dc "$p-two-moment-charge-from$c-x05-checkpoint-$d.json.gz" > "/tmp/remap-archive-from$c.json"
done
python -W error early_pole_endpoint_charge_grid_donor_terminal_receipt.py results \
  "$p-endpoint-crossed-charge-terminal-complete-$d.json" \
  --records "$p-two-moment-charge-from2-x05-receipt-$d.json" "$p-two-moment-charge-from4-x05-receipt-$d.json" \
  --payloads /tmp/remap-archive-from{2,4}.json
```
