# A positive local pressure discrepancy can become negative after propagation

Use the exact evolved G9 x=.04 start. For the first small interval of width log(1.25)/1024, compare one update with two half steps using the same finite-capacity step. The resulting two states exactly match the earlier embedded-step probe. From that instant, both states take the same 512-step uniform log-time continuation to x=.05. The first interval is the only intentionally different operation. Every pressure here is evaluated at the same native x=.05 measurement mass, including the near-start comparison.

| Two half steps minus one update, relative to the one-update arm | Immediately after the first interval | At x=.05 after identical continuation |
|---|---:|---:|
| Total count | +3.4156960e-7 | +5.2108777e-9 |
| Total pressure | +6.3262627e-10 | −2.4517194e-10 |
| Warm pressure | +1.7510734e-10 | −2.5017325e-10 |

Cold cells contribute 72.36% of the first absolute cellwise pressure discrepancy. At the endpoint the warm share is 96.78% and the signed warm contribution is −2.4977e-10. Number differences remain cold-dominated. The sign change gives a direct, controlled counterexample to using an unpropagated local pressure difference as a terminal pressure estimate. It does not decompose all errors in the full shifted grid, since only one step was changed.

Both completed 512-step payloads retain the exact start-dependent signature and the same original positivity, pair-number and continuity gates. They replay the stored result byte-for-byte without taking new continuation steps.

```bash
cd kinetics
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
d=2026-11-28; p=results/early-pole-fine-charge2-scalar4-gauss9warm; w=results/early-pole-warm
tmp=$(mktemp -d)
gzip -dc "$w-gauss9-propagated-step-one-checkpoint-$d.json.gz" > "$tmp/one.json"
gzip -dc "$w-gauss9-propagated-step-two-checkpoint-$d.json.gz" > "$tmp/two.json"
python -W error early_pole_endpoint_warm_propagated_embedded_defect.py "$p-x04-state-$d.json.gz" \
 "$w-checkpoint-clock-comparison-x05-$d.json" "$w-embedded-step-probe-$d.json" "$tmp/one.json" "$tmp/two.json"
rm -r "$tmp"
```

Finite comparison on one numerical G9 chart only. No continuum accuracy or physical cosmological inference.
