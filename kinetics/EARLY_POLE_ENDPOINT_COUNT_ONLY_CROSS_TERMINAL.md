# Swapping only the initial cold count reproduces nearly all of the scalar donor effect

The full D2/D4 crossed terminal histories separated the charged donor from the inherited scalar donor. Now rerun both .04→.05 controls after keeping **each charged chart's native scalar shape** throughout the initial state, except for one multiplicative rescale of the cold q≤1e−5 group to give it the other donor's cold count. The outside scalar cells and all charged occupations stay native at the restart. Continue the same 2048-step clock and recheck complete advanced occupations rather than extrapolating a linear gap estimate.

| Outcome | Path-averaged count-only replacement | Path-averaged full scalar replacement | Count-only share |
|---|---:|---:|---:|
| Terminal cold count | −8.84993484e−5 | −8.90812241e−5 | 99.347% |
| Terminal all-scalar count | −8.84991079e−5 | −8.90796074e−5 | 99.348% |
| Terminal scalar Hartree pole | −1.49445634e−6 | −1.50427015e−6 | 99.348% |

The individual terminal cold-count differences between the count-only and full swaps are +6.461e−7 on native D2 charges and −5.177e−7 on native D4 charges. Do not interpret a 99.35% share of the **scalar donor effect** as 99.35% of the native D2→D4 discrepancy; the separate charged donor was larger. The small residual combines shape and outside-cold differences and does not certify physical insensitivity to all scalar shapes.

```bash
cd kinetics
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
d=2026-11-28
p=results/early-pole
h="$p-endpoint-charge-cold-count-history-budget-$d.json"
for c in 2 4; do
  s=$((6-c))
  python -W error early_pole_endpoint_count_only_cross_terminal.py results "$h" \
    "$p-crossed-charge-$c-scalar-$s-x05-receipt-$d.json" "$c" "/tmp/count-only-replay-C$c.drive.json"
  gzip -dc "$p-count-only-cross-charge-$c-x05-checkpoint-$d.json.gz" > "/tmp/count-only-C$c.drive.json"
done
python -W error early_pole_endpoint_count_only_cross_terminal_receipt.py results "$h" \
  "$p-endpoint-crossed-charge-terminal-complete-$d.json" \
  --full-cross "$p-crossed-charge-2-scalar-4-x05-receipt-$d.json" "$p-crossed-charge-4-scalar-2-x05-receipt-$d.json" \
  --records "$p-count-only-cross-charge-2-x05-receipt-$d.json" "$p-count-only-cross-charge-4-x05-receipt-$d.json" \
  --payloads /tmp/count-only-C{2,4}.drive.json
```
