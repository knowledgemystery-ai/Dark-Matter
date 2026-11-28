# In crossed .04→.05 trajectories the charged donor dominates the cold-count contrast

At x=.04 keep the native charged populations and charge grid from one of the D2/D4 independently evolved branches, and take only the scalar occupation vector from either branch. The scalar grid is exactly common, so no remap is involved. The mixed starts are conditional controls, not pre-existing physical histories. Evolve both mixed starts through the original 2048-step terminal clock; the original native endpoints supply the diagonal.

| Charged chart and occupations | D2 initial scalar | D4 initial scalar |
|---|---:|---:|
| D2 | .01738683682 | .01728994990 |
| D4 | .01697933293 | .01689805741 |

The table gives cold scalar counts at x=.05 for q≤1e−5. The native D2→D4 diagonal change is −.000488779413. The charged-donor replacement averaged along the two paths around the square contributes −.000399698189 (**81.77%** of the diagonal change), while the scalar-donor replacement contributes −.000089081224 (**18.23%**). Their sum is the diagonal. The interaction is +.000015611395, measured separately; it is already reflected in the path averages and must not be added to them. Both single replacements decrease the endpoint on either path. The relative allocation for the scalar pole shift is 81.78% and 18.22%.

The maximum absolute pair-count residual is 3.30e−17, and the maximum continuity residual in either mixed replay is 8.29e−12. The original native trajectories are independently replayed from their same .04 inputs with the same step rule; both reproduce every endpoint occupation bitwise. Exact final occupation checkpoints are retained for both mixed arms and the reducer independently reconstructs each pole and count.

This only separates the inherited scalar state from the combined charged **grid and occupations** under this one late restart. It neither disentangles the grid from the inherited charged occupations nor assigns a unique cause to the full earlier evolution. It does not assess a physical production rate or the separate cold photon-assisted band convergence.

## Reproduce the exact reduction

```bash
cd kinetics
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
d=2026-11-28
p=results/early-pole
for c in 2 4; do
  s=$((6-c))
  gzip -dc "$p-crossed-charge-$c-scalar-$s-x05-checkpoint-$d.json.gz" > "/tmp/cross-final-C$c-S$s.json"
  python -W error early_pole_endpoint_native_terminal_diagonal_audit.py results "$c"
done
python -W error early_pole_endpoint_crossed_charge_terminal_receipt.py results \
  "$p-endpoint-charge-cold-count-history-budget-$d.json" \
  "$p-endpoint-charge-cold-count-gap-ladder-$d.json" \
  --records "$p-crossed-charge-2-scalar-4-x05-receipt-$d.json" "$p-crossed-charge-4-scalar-2-x05-receipt-$d.json" \
  --final-checkpoints /tmp/cross-final-C{2-S4,4-S2}.json
```
