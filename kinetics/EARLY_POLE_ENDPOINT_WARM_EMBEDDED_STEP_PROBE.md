# Local embedded defects are not the global time-grid error

Take the exact independently accepted G9 and mixed300 states at x=.02,.04,.05. From each hold point evolve over a single small log-time interval using one, two or four equal updates and the same capacity/Pauli step. Repeat for a second interval four times wider. Each pressure measurement keeps the common native mass belonging to that **starting** hold point; no physical expansion of the measurement mass is inferred. The final x=.05 pulse is only a local numerical probe, not a new physical endpoint.

For the smaller interval, log(1.25)/1024, G9 one-to-two relative count changes are +3.4241e-9, +3.4187e-7 and +7.2254e-7 at the three states; pressure changes are +1.7356e-11, +7.6545e-10 and +1.4810e-9. Mixed300 gives nearly the same count and pressure changes. The relative one-to-two pressure change rises from .04 to .05 along with count; the walk-forward shifted-grid benchmark nevertheless improved count and worsened terminal pressure. Local endpoint trends alone therefore do not establish a pressure-accuracy policy.

At x=.04 and x=.05, cold cells contribute about 77% and 74% of the absolute local pressure change, whereas around 80% of the **accumulated** shifted-versus-reference pressure difference on the earlier archived .02–.04 interval came from warm cells. The local count change remains almost entirely cold at all three states. These are differences between observables measured at different parts of the procedure, not an inferred causal decomposition of the global error. Signed cell sums close the direct local measurement and paired-particle changes retain their prior conservation and positivity gates.

```bash
cd kinetics
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
d=2026-11-28; w=results/early-pole-warm
python -W error early_pole_endpoint_warm_embedded_step_probe.py results "$w-single-g9-upgrade-$d.json" \
 "$w-checkpoint-clock-comparison-x02-$d.json" "$w-checkpoint-clock-comparison-x04-$d.json" \
 "$w-checkpoint-clock-comparison-x05-$d.json"
```

This is a finite-context numerical diagnostic, not a global error estimator, asymptotic order proof or physical matter prediction.
