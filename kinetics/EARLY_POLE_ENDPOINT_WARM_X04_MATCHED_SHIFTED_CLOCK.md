# Matched shifted 1024-step G9 and mixed300 clocks

Use the same predeclared nonuniform log-time map in both charts, with exponent 1.2, unchanged reaction step and identical endpoint times x=.02 and .04. No mixed300 coefficient was refitted. Each calculation starts with its own exact accepted checkpoint and is compared with its own uniform 4096-step result at the common native pressure mass.

| Discrepancy relative to uniform 4096 | Uniform 1024 | Shifted 1024 | Shifted/uniform magnitude |
|---|---:|---:|---:|
| G9 count | −1.59383e-4 | −6.63544e-5 | 0.4163 |
| mixed300 count | −1.59383e-4 | −6.63544e-5 | 0.4163 |
| G9 pressure | +3.83234e-7 | +3.41309e-7 | 0.8906 |
| mixed300 pressure | +3.75170e-7 | +3.36864e-7 | 0.8979 |
| mixed300/G9 count contrast | +9.63272e-13 | +9.35980e-14 | 0.0972 |
| mixed300/G9 pressure contrast | −8.06387e-9 | −4.44531e-9 | 0.5513 |

For the individual observables, the shifted 1024-step result does not beat uniform 2048. The pressure contrast does not beat uniform 2048 either. The smaller count contrast is a cancellation and not a numerical error certificate. The paired comparison computes small differences by the stable ratio identity before subtracting any nearly equal overall ratios. In each chart practically all absolute shifted count discrepancy comes from cold cells, while about 80% of absolute pressure discrepancy comes from warm cells. The leading cold cell is q in [3.306e-6,3.625e-6] and the leading pressure cell is [0.62843,0.90829]. An independent replay of the three uniform clocks confirms that the count contrast's change is at most 2.072e-8 of the larger individual shift at **all eight sampled outputs**, not merely the last.

These are actual same-start numerical checks, not estimates of continuum error, a physical abundance or a transfer observable.

```bash
cd kinetics
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
d=2026-11-28; p=results/early-pole-fine-charge2-scalar4; w=results/early-pole-warm
tmp=$(mktemp); trap 'rm -f "$tmp"' EXIT
gzip -dc "$w-mixed300-x04-shifted1024-checkpoint-$d.json.gz" > "$tmp"
python -W error early_pole_endpoint_warm_mixed300_prescribed_clock.py \
 "$p-mixed300warm-x02-state-$d.json.gz" "$p-mixed300warm-x04-state-$d.json.gz" \
 "$w-mixed300-x04-third-clock-$d.json" "$w-single-g9-upgrade-$d.json" "$tmp"
python -W error early_pole_endpoint_warm_x04_paired_shifted_clock.py \
 "$w-gauss9-x04-shifted1024-$d.json" "$w-mixed300-x04-shifted1024-$d.json" \
 "$w-gauss9-x04-third-clock-$d.json" "$w-mixed300-x04-third-clock-$d.json" \
 "$w-x04-third-clock-paired-contrast-$d.json"
python -W error early_pole_endpoint_warm_x04_paired_clock_count_history.py \
 "$p-gauss9warm-x04-$d.json.gz" "$p-mixed300warm-x04-$d.json.gz" \
 "$w-gauss9-x04-third-clock-$d.json" "$w-mixed300-x04-third-clock-$d.json" \
 "$w-checkpoint-clock-comparison-x04-$d.json"
python -W error early_pole_endpoint_warm_shifted_clock_partition.py \
 "$p-mixed300warm-x04-state-$d.json.gz" "$w-mixed300-x04-shifted1024-checkpoint-$d.json.gz" \
 "$w-mixed300-x04-4096-checkpoint-$d.json.gz" "$w-mixed300-x04-shifted1024-$d.json" \
 "$w-mixed300-x04-third-clock-$d.json" --control "$w-single-g9-upgrade-$d.json"
```

For the G9 partition, substitute the corresponding G9 inputs and omit the control option.
