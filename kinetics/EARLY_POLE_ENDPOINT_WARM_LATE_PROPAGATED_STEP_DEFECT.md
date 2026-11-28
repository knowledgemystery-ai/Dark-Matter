# Repeating the propagated pressure sign reversal at a held interior start

Repeat the controlled one-step versus two-half-step substitution, now starting from the actual G9 uniform4096 checkpoint at prefix 2560, x=0.04598633. The checkpoint's full payload hash, signature and ledger prefixes match its completed uniform reference. The pulse width log(1.25)/1024 and the native x=.05 pressure mass are exactly as in the earlier first-step experiment. After the pulse each arm uses the same 192-step continuation to x=.05. A continuation step is 0.998371 times the width used in the earlier 512-step continuation, so future evolution is tested at comparable time resolution.

| Relative two-minus-one difference | At the end of the pulse | At x=.05 |
|---|---:|---:|
| Count | +2.91963e-7 | +2.41669e-8 |
| Pressure | +6.03403e-10 | −2.70111e-10 |
| Warm pressure | +1.68166e-10 | −3.06192e-10 |

The same pressure sign reversal therefore occurs for this later, independent pulse. Warm cells make up 27.83% of the absolute pressure discrepancy immediately after the pulse and 90.03% after propagation. The earlier first pulse gave +6.32626e-10 then −2.45172e-10 respectively. This does not establish that every insertion reverses sign or decompose the error of a whole time grid.

Both actual paired experiments meet the original number, positivity and continuity gates. The completed two-arm checkpoints are fingerprinted and replay their measurements without additional continuation steps.

```bash
cd kinetics
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
d=2026-11-28; p=results/early-pole-fine-charge2-scalar4-gauss9warm; w=results/early-pole-warm
tmp=$(mktemp -d)
for a in one two; do
 gzip -dc "$w-gauss9-late-propagated-step-$a-checkpoint-$d.json.gz" > "$tmp/$a.json"
done
python -W error early_pole_endpoint_warm_late_propagated_defect.py "$p-x05-state-$d.json.gz" \
 "$w-gauss9-x05-uniform4096-progress-2560-$d.json.gz" "$w-gauss9-x05-uniform4096-checkpoint-$d.json.gz" \
 "$w-gauss9-x05-uniform4096-$d.json" "$w-gauss9-propagated-step-defect-$d.json" "$tmp/one.json" "$tmp/two.json"
rm -r "$tmp"
```
