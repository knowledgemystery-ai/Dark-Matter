# Embedded local defects on actual interior x05 reference states

Use the exact x=.04 and .05 endpoints plus five in-flight checkpoints saved at common predetermined prefix steps 1536, 2048, 2560, 3072 and 3584 of each 4096-step run. A continued checkpoint is accepted only if its full payload hash, signature, time, and all recorded work, caps and ledger prefixes match the completed reference. No intermediate occupation is interpolated.

At each held state independently compare one, two and four steps over the same two small intervals. Pressure measurements for **all** sample locations are frozen at the same native x=.05 mass. The smaller interval, log(1.25)/1024, gives the following G9 two-minus-one differences relative to the held totals:

| Exact uniform4096 prefix | x | Count | Pressure | Cold share of local absolute pressure change |
|---:|---:|---:|---:|---:|
| 0 | 0.04000000 | +3.418709e-07 | +6.330326e-10 | 72.360% |
| 1536 | 0.04349119 | +3.541103e-07 | +7.440916e-10 | 72.088% |
| 2048 | 0.04472136 | +4.264246e-07 | +8.511580e-10 | 72.299% |
| 2560 | 0.04598633 | +2.921448e-07 | +6.037864e-10 | 72.173% |
| 3072 | 0.04728708 | +8.441537e-07 | +1.650285e-09 | 75.148% |
| 3584 | 0.04862462 | +6.457114e-07 | +1.320596e-09 | 73.377% |
| 4096 | 0.05000000 | +7.230696e-07 | +1.481251e-09 | 74.046% |

The mixed300 local pressure differences follow the same pattern; their largest absolute difference from the G9 relative pulse, across these states, is 9.237e-12. Compare these with the separately published accumulated clock discrepancies: accumulated pressure differences need not share the local cell attribution. This profile measures actual local truncation contrasts, not their propagated effect on a terminal observable, a continuum bound, or a physical cosmological prediction.

```bash
cd kinetics
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
d=2026-11-28; p=results/early-pole-fine-charge2-scalar4; w=results/early-pole-warm
for a in gauss9 mixed300; do
 python -W error early_pole_endpoint_warm_x05_interior_clock_probe.py "$p-${a}warm-x04-state-$d.json.gz" "$p-${a}warm-x05-state-$d.json.gz"   "$w-$a-x05-uniform4096-checkpoint-$d.json.gz" "$w-$a-x05-uniform4096-$d.json" "$w-checkpoint-clock-comparison-x05-$d.json"   "$w-single-g9-upgrade-$d.json" "$w-$a-x05-uniform4096-progress-"{1536,2048,2560,3072,3584}"-$d.json.gz"
done
```

The observed profile is not monotone between the old endpoints. In G9 the relative pressure pulse drops from 8.512e-10 at prefix 2048 to 6.038e-10 at 2560, then jumps to 1.650e-9 at 3072. Thus interpolating solely from the beginning and end would miss an interior dip and peak even before asking how these defects propagate.
