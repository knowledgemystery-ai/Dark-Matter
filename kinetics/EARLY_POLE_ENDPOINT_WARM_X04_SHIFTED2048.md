# Actual matched 2048-step shifted refinement

Repeat the same predeclared nonuniform log-time map with exponent 1.2 at 2048 steps for G9 and mixed300, starting .02 and ending .04. Each chart uses its own exact archived starting occupation, unchanged spatial grid and gain/loss update; the refined shift is not extrapolated from the 1024 calculation. All measurements use the same native reference pressure mass.

| Observable | Uniform 2048 relative to uniform 4096 | Shifted 2048 relative to uniform 4096 | Shifted/uniform discrepancy magnitude |
|---|---:|---:|---:|
| G9 count | -3.1715696e-05 | -8.8238483e-06 | 0.27822 |
| mixed300 count | -3.1715695e-05 | -8.8238483e-06 | 0.27822 |
| G9 pressure | +8.1709524e-08 | +7.2423461e-08 | 0.88635 |
| mixed300 pressure | +7.9463592e-08 | +7.1243039e-08 | 0.89655 |
| mixed300/G9 count contrast | +5.8455096e-13 | +6.6058853e-14 | 0.11301 |
| mixed300/G9 pressure contrast | -2.2458420e-09 | -1.1803751e-09 | 0.52558 |

The contrast row is a change relative to the accepted 4096 contrast, not an individual error. Exact completed step checkpoints retain trajectory checks and original conservation gates in both charts. This is a comparison to one available fine clock, not an asymptotic order measurement, continuum estimate, or physical prediction.

```bash
cd kinetics
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
d=2026-11-28; p=results/early-pole-fine-charge2-scalar4; w=results/early-pole-warm
for a in gauss9 mixed300; do
 tmp=$(mktemp)
 gzip -dc "$w-$a-x04-shifted2048-checkpoint-$d.json.gz" > "$tmp"
 args=()
 if [ "$a" = mixed300 ]; then args=(--control "$w-single-g9-upgrade-$d.json"); fi
 python -W error early_pole_endpoint_warm_x04_shifted2048.py "$p-${a}warm-x02-state-$d.json.gz" "$p-${a}warm-x04-state-$d.json.gz" "$w-$a-x04-third-clock-$d.json" "$w-$a-x04-shifted1024-$d.json" "$tmp" "${args[@]}"
 rm -f "$tmp"
done
python -W error early_pole_endpoint_warm_x04_shifted_refinement_pair.py "$w-gauss9-x04-shifted2048-$d.json" "$w-mixed300-x04-shifted2048-$d.json" "$w-x04-paired-shifted-clock-$d.json" "$w-x04-third-clock-paired-contrast-$d.json"
```
