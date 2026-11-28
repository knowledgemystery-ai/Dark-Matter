# One last held pressure pulse close to x=.05

Use the genuine G9 uniform4096 prefix checkpoint at 3584 (x=0.04862462), with exact payload and full-ledger prefix binding to the accepted final state. Hold the same x05 measurement mass and pulse width log(1.25)/1024 as in the earlier experiments. Both arms then take the same 64-step continuation; the width differs by less than 0.7% from the previous first-pulse experiment's continuation width.

| Two half steps minus one, relative difference | Just after the pulse | At x=.05 |
|---|---:|---:|
| Count | +6.4523894e-7 | +4.5811360e-8 |
| Pressure | +1.3197607e-9 | −1.0757944e-10 |
| Warm pressure | +3.5191792e-10 | −1.6034107e-10 |

The observed pressure sign reversal is reproduced even from this nearer endpoint, although its magnitude changes; the warm share of the absolute pressure difference shifts from about 26.6% to 72.4%. Along with the two earlier insertions this makes three actual controlled counterexamples to using the local pressure sign as the terminal sign. It does not prove this sign behavior at other insertion times or the sign of a summed grid discrepancy.

```bash
cd kinetics
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
d=2026-11-28; p=results/early-pole-fine-charge2-scalar4-gauss9warm; w=results/early-pole-warm
tmp=$(mktemp -d)
for a in one two; do gzip -dc "$w-gauss9-near-end-propagated-step-$a-checkpoint-$d.json.gz" > "$tmp/$a.json"; done
python -W error early_pole_endpoint_warm_near_end_propagated_defect.py "$p-x05-state-$d.json.gz" \
 "$w-gauss9-x05-uniform4096-progress-3584-$d.json.gz" "$w-gauss9-x05-uniform4096-checkpoint-$d.json.gz" \
 "$w-gauss9-x05-uniform4096-$d.json" "$w-gauss9-propagated-step-defect-$d.json" "$tmp/one.json" "$tmp/two.json"
rm -r "$tmp"
```
