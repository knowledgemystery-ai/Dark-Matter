# Does native pressure determine the common-phase frozen streaming sign?

Retain three quantities independently in each refined warm shell of the exact x=.05 G8, G9 and mixed300 states: number, mean excitation energy, and pressure at the common reference mass. Keep the cold convention unchanged. All choices use the same phase from the archived G8 nodal case for each trial interval.

With these three constraints the result for A=2 is again ordered G9 above mixed300 above G8; in particular the mixed300−G9 outer interval is approximately [−4.913e-11,−4.026e-11], where the previous number-and-energy-only comparison allowed both signs. For A=10, 100 and 1000 every pair still has positive measures giving **both signs**. For A=10, G9−G8 constructive contrasts include −1.4606e-10 and +1.1545e-9. The exact positive atom supports and weights for those two comparisons are included, and an independent verifier checks them without solving an optimization problem. These are not statements about physical scattering or a calibrated cosmic scale.

At fixed number and energy, imposing pressure q²/E = E−M²/E is equivalent to fixing the inverse-energy moment. Write E(t)=E0+h t in a normalized shell coordinate t; subtract the constant and linear parts of 1/E to use the well-scaled feature t(1−t) E0/E(t). Linear optimization on a grid containing the original atoms then gives three-point positive witnesses. A dual affine combination of these three features, corrected for its measured mesh residual and a controlled interpolation remainder, supplies a continuous outer bound. The nonlinear feature has second derivative −2 E0²(E0+h)/E(t)³, enabling the same direct interpolation-error estimate as for the velocity kernel. The final outer bound intersects the previous number-and-energy-only enclosure because the feasible set is a subset. It uses ordinary floating-point arithmetic and numerically integrated fixed cold cells, not directed rounding. Independent witness checks verify number, mean energy and the conditioned pressure feature rather than relying on the optimizer tolerance.

```bash
cd kinetics
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
d=2026-11-28; p=results/early-pole-fine-charge2-scalar4; w=results/early-pole-warm
python -W error -m unittest -v test_early_pole_endpoint_warm_streaming_pressure_constraint
python -W error early_pole_endpoint_warm_streaming_pressure_constraint.py "$p-x05-state-$d.json.gz" \
 "$w-mixed300-actual-prefix-x05-$d.json" "$w-single-g9-upgrade-$d.json" \
 "$w-three-chart-radiation-streaming-$d.json" "$w-three-chart-streaming-hull-$d.json" \
 "$p"-{gauss8warm,gauss9warm,mixed300warm}-x05-state-$d.json.gz
python -W error early_pole_endpoint_warm_streaming_pressure_witness.py "$p-x05-state-$d.json.gz" \
 "$w-three-chart-pressure-constrained-streaming-$d.json" "$p"-gauss{8,9}warm-x05-state-$d.json.gz
```
