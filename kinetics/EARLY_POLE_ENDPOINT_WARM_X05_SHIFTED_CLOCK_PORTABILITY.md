# Out-of-interval test of the fixed shifted clock

Use the identical log-time schedule exponent 1.2 from x=.02→.04, without retuning, on the independent shorter x=.04→.05 interval. Both analyses use exact native archived starts and the same endpoint reference mass. At fixed 1024 steps compare with the original uniform 1024, using the independently evolved uniform 2048 as the available reference. These are new trajectories, not an extrapolation of the longer interval.

| Observable | Uniform 1024 minus reference, fraction | Shifted 1024 minus reference, fraction | Shifted/uniform discrepancy magnitude |
|---|---:|---:|---:|
| G9 count | -3.1800957e-05 | -2.6204753e-05 | 0.824024 |
| mixed300 count | -3.1800956e-05 | -2.6204753e-05 | 0.824024 |
| G9 pressure | +2.5244576e-07 | +5.2261084e-07 | 2.070191 |
| mixed300 pressure | +2.4896116e-07 | +5.1284035e-07 | 2.059921 |
| Paired count contrast | +1.9851419e-13 | -9.3150153e-14 | 0.469237 |
| Paired pressure contrast | -3.4844536e-09 | -9.7700778e-09 | 2.803905 |

In particular a step prescription must be assessed separately for individual counts, individual pressure, and the residuals of the paired contrast. The smaller paired count does not bound the individual. The observed changes are local reference discrepancies, not a forecast for other starts, a continuum solution, a transfer estimate or physical dark-matter abundance.

```bash
cd kinetics
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
d=2026-11-28; p=results/early-pole-fine-charge2-scalar4; w=results/early-pole-warm
for a in gauss9 mixed300; do
 tmp=$(mktemp)
 gzip -dc "$w-$a-x05-shifted1024-checkpoint-$d.json.gz" > "$tmp"
 args=(); if [ "$a" = mixed300 ]; then args=(--control "$w-single-g9-upgrade-$d.json"); fi
 python -W error early_pole_endpoint_warm_shifted_x05.py "$p-${a}warm-x04-state-$d.json.gz" "$p-${a}warm-x05-state-$d.json.gz" "$w-checkpoint-clock-comparison-x05-$d.json" "$tmp" "${args[@]}"
 rm -f "$tmp"
done
python -W error early_pole_endpoint_warm_shifted_x05_pair.py "$w-gauss9-x05-shifted1024-$d.json" "$w-mixed300-x05-shifted1024-$d.json" "$w-checkpoint-clock-comparison-x05-$d.json"
```

The shorter shifted count discrepancy improves only about 17.6%, unlike the roughly 58% gain on the longer held restart with the finer reference. The individual pressure discrepancy approximately doubles, and the pressure contrast grows by about 2.80-fold. Thus this particular prescription is **not** a generally better count-and-pressure time grid. There is no retuned retry in this experiment.
