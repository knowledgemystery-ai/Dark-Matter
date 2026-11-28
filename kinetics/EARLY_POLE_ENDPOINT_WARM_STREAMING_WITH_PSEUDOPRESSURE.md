# Retaining native pseudo-pressure improves but does not fix frozen streaming

Reapply the actual common-phase G8/G9/mixed300 x=.05 comparison, with the original cold convention. Independently hold number, mean excitation energy, native pressure, and native pseudo-pressure in each of the 22 warm shells. All charts and shapes share the same fixed G8 nodal phase at each of the prescribed expansion intervals.

At A=10 the four constraints now ensure G9 and mixed300 both give higher loss than G8. For G9−G8 the outer range is about [+1.486e-10,+8.599e-10]; three-moment sign witnesses had previously spanned both signs. However, mixed300−G9 still has explicit witnesses at roughly −3.393e-10 and +2.446e-10 at A=10. For A=100 and 1000, both signs are still attainable in **every** comparison. The A=100 G9/G8 constructive extremes (about −8.002e-9,+9.035e-9) are also serialized as positive atom supports and verified separately without optimizing.

At fixed number and energy, ordinary pressure and pseudo-pressure are equivalent to retaining the inverse energy and inverse cube moments respectively. The algebra is q²/E=E−M²/E and q⁴/E³=E−2M²/E+M⁴/E³. In normalized shell energy t with E=E0(1+r t), combine inverse cube with the existing constant, linear and inverse features. The stable fourth feature is

```text
K(t) = t^3(4+3rt)/(1+rt)^3
     = [6/(1+rt)-5+3rt-(1+rt)^(-3)]/r^3.
```

The rational feature has second derivative 12t(2+rt)/(1+rt)^5, bounded by 12(2+r). This supplies the interpolation remainder in the dual enclosure. The optimizer chooses at most four positive atoms per cell. The output intersects the coarser available three-moment bounds, and tests check algebraic equivalence and independent positive null-direction perturbations in both pressure moments. Bounds use ordinary floating point and numerical quadrature for the unchanged cold cells.

These are formal frozen transport probes, not a physical transfer function or a prediction about an observed dark matter spectrum.

```bash
cd kinetics
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
d=2026-11-28; p=results/early-pole-fine-charge2-scalar4; w=results/early-pole-warm
python -W error -m unittest -v test_early_pole_endpoint_warm_streaming_pseudopressure_constraint
python -W error early_pole_endpoint_warm_streaming_pseudopressure_constraint.py "$p-x05-state-$d.json.gz" \
 "$w-mixed300-actual-prefix-x05-$d.json" "$w-single-g9-upgrade-$d.json" \
 "$w-three-chart-radiation-streaming-$d.json" "$w-three-chart-pressure-constrained-streaming-$d.json" \
 "$p"-{gauss8warm,gauss9warm,mixed300warm}-x05-state-$d.json.gz
python -W error early_pole_endpoint_warm_streaming_pseudopressure_witness.py "$p-x05-state-$d.json.gz" \
 "$w-three-chart-pseudopressure-constrained-streaming-$d.json" "$p"-gauss{8,9}warm-x05-state-$d.json.gz
```

## Effect on the exact G9/G8 shape budget

The minimum number of shared warm cells whose actual positive shape must be pinned moves from 0,0,1,10,14 under the three constraints to **0,0,0,9,14** under the four constraints at A=1,2,10,100,1000. The corresponding constructive insufficiency and outer sufficient cardinalities meet at all five points on the 1536-subinterval mesh. In particular the shortest ambiguous reference no longer needs an extra exact shape after holding pseudo-pressure; the longest needs just as many shapes as before. This measures a worst-case finite-moment budget at the selected phases, not the cost of evolving the cells.

```bash
python -W error early_pole_endpoint_warm_streaming_pseudopressure_pin_budget.py "$p-x05-state-$d.json.gz" \
 "$w-three-chart-pseudopressure-constrained-streaming-$d.json" "$p"-gauss{8,9}warm-x05-state-$d.json.gz
```
