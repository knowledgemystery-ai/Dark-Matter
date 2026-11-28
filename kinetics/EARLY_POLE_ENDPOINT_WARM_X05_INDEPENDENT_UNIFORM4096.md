# Independent 4096-step check of short-window shifted-clock portability

Compute a new genuine uniform 4096-step trajectory for each selected warm chart from its own exact x=.04 state. Compare all trajectories at x=.05 and at the same native measurement mass. This tests the earlier comparison against uniform2048 without retuning, extrapolating a new solution, or comparing unequal start states.

| Observable | Uniform 1024 vs uniform 4096 | Shifted 1024 vs uniform 4096 | Uniform 2048 vs uniform 4096 | Shifted/uniform 1024 discrepancy |
|---|---:|---:|---:|---:|
| G9 count | -3.9568979e-05 | -3.3972819e-05 | -7.7682697e-06 | 0.85857 |
| mixed300 count | -3.9568979e-05 | -3.3972819e-05 | -7.7682699e-06 | 0.85857 |
| G9 pressure | +3.1520967e-07 | +5.8537477e-07 | +6.2763897e-08 | 1.85710 |
| mixed300 pressure | +3.1094648e-07 | +5.7482569e-07 | +6.1985308e-08 | 1.84863 |
| Paired count contrast | +5.5957455e-14 | -2.3570836e-13 | -1.4266477e-13 | 4.21228 |
| Paired pressure contrast | -4.2630095e-09 | -1.0548634e-08 | -7.7855566e-10 | 2.47446 |

The individual count discrepancy continues to improve with the shift; the pressure discrepancy continues to worsen relative to now-independent uniform4096. The paired entries are signed differences from the 4096 ratio, not individual numerical errors. The earlier signed ratios against the 2048 reference reproduce in the same audit. Exact completed checkpoints use the same conservation and positivity gates. These are finite-reference numerical comparisons, not continuum accuracy, a physical abundance, or cosmological transfer.

```bash
cd kinetics
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
d=2026-11-28; p=results/early-pole-fine-charge2-scalar4; w=results/early-pole-warm
for a in gauss9 mixed300; do
 tmp=$(mktemp)
 gzip -dc "$w-$a-x05-uniform4096-checkpoint-$d.json.gz" > "$tmp"
 args=()
 if [ "$a" = mixed300 ]; then args=(--control "$w-single-g9-upgrade-$d.json"); fi
 python -W error early_pole_endpoint_warm_x05_uniform4096.py "$p-${a}warm-x04-state-$d.json.gz" "$p-${a}warm-x05-state-$d.json.gz" "$w-checkpoint-clock-comparison-x05-$d.json" "$w-$a-x05-shifted1024-$d.json" "$tmp" "${args[@]}"
 rm -f "$tmp"
done
python -W error early_pole_endpoint_warm_x05_uniform4096_pair.py "$w-gauss9-x05-uniform4096-$d.json" "$w-mixed300-x05-uniform4096-$d.json" "$w-checkpoint-clock-comparison-x05-$d.json" "$w-x05-shifted-portability-$d.json"
```

The pressure deterioration survives the finer reference, though its factor changes: about 1.86 and 1.85 in the two charts, compared with about 2.07 and 2.06 against 2048. The paired count conclusion **does change**: its shifted/unshifted discrepancy factor goes from about 0.47 against 2048 to about 4.21 against 4096, with a sign reversal. These count contrast differences remain extremely small (order 1e-13), not estimates of the individual count error.

Against the 4096 reference, the cold cells contain essentially all of the individual count discrepancy. For G9, warm cells contribute 86.2% of the absolute uniform1024 pressure discrepancy and 92.7% of the shifted pressure discrepancy, with cold and warm signed contributions opposing. The analogous mixed300 fractions are 86.1% and 92.6%. In all four pressure comparisons the largest cellwise contributor is q in [0.90829,1.31278]. These statements are additive measurements of the actual endpoint differences, not a causal decomposition of the method error.

```bash
# Same directory and variables as above.
for a in gauss9 mixed300; do
 args=()
 if [ "$a" = mixed300 ]; then args=(--control "$w-single-g9-upgrade-$d.json"); fi
 python -W error early_pole_endpoint_warm_x05_clock_cell_attribution.py "$p-${a}warm-x05-state-$d.json.gz" \
  "$w-$a-x05-1024-checkpoint-$d.json.gz" "$w-$a-x05-shifted1024-checkpoint-$d.json.gz" \
  "$w-$a-x05-uniform4096-checkpoint-$d.json.gz" "$w-checkpoint-clock-comparison-x05-$d.json" \
  "$w-$a-x05-shifted1024-$d.json" "$w-$a-x05-uniform4096-$d.json" "${args[@]}"
done
```
