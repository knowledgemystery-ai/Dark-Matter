# One shifted 1024-step clock reduces the same-start G9 count discrepancy

Continue the exact actual x=.02 G9 state to x=.04 with the same geometry, finite-capacity reaction update, and 1024 steps. Only the log-time spacing changes. For u_i=i/1024, use t_i/log(2)=log(1+u_i (exp(1.2)-1))/1.2. The coefficient was fixed before running this shifted candidate, based on the steadily increasing count discrepancy in the archived uniform-grid segment. The method spends more updates late; the actual interval widths range from 0.5825 to 1.9312 times the 1024-step uniform width. This is one prescribed comparison, not an error controller fitted to the output.

Relative to the same-start **uniform 4096-step** reference, the absolute scalar count discrepancy decreases from 1.59383e-4 (uniform 1024) to 6.63544e-5 (shifted 1024), about 58.37% smaller at the same step count. The corresponding pressure discrepancies are 3.83234e-7 and 3.41309e-7, about 10.94% smaller. The shifted 1024 result does not beat uniform 2048 in either quantity; these are reference discrepancies, not errors against an exact continuum solution.

The ordinary continuity, paired-number and Pauli gates are unchanged. The exact completed step checkpoint binds the G9 geometry, restart, kernel and full shifted time vector, and its no-step replay reproduces the output byte-for-byte. This does not establish improvement on other clocks, warm orders or the full interval from startup.

```bash
cd kinetics
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
d=2026-11-28; p=results/early-pole-fine-charge2-scalar4-gauss9warm; w=results/early-pole-warm
tmp=$(mktemp); trap 'rm -f "$tmp"' EXIT
gzip -dc "$w-gauss9-x04-shifted1024-checkpoint-$d.json.gz" > "$tmp"
python -W error early_pole_endpoint_warm_g9_prescribed_clock.py "$p-x02-state-$d.json.gz" \
 "$p-x04-state-$d.json.gz" "$w-gauss9-x04-third-clock-$d.json" "$tmp"
```

For a separate direct replay instead of reopening the completed checkpoint, start with an absent temporary checkpoint path.
