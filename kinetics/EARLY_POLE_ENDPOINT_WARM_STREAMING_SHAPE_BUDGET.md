# Which warm shapes are needed to preserve the common-phase sign?

In the frozen G9/G8 comparison at the exact x=.05 endpoint, cold cells remain fixed. In each of the 22 warm cells the archived number, energy and native pressure remain fixed. As a counterfactual, choose common momentum cells and hold their entire **actual nodal shape** fixed in both charts. Ask how many selected cells are needed to guarantee the sign of the archived nodal contrast at the fixed G8 phase, allowing all admissible variations in the other cells.

| Prescribed A | Nodal sign (G9−G8) | Fewest cells that must be pinned |
|---:|:---:|---:|
| 1 | + | 0 |
| 2 | + | 0 |
| 10 | + | 1 |
| 100 | − | 10 |
| 1000 | + | 14 |

At A=10, the single sufficient cell is q in [0.01579697, 0.02283191], with a guaranteed positive margin at least 1.2588e-10 in absolute loss. This is different from the cell that contributed most of the time-step pressure difference; the observable and the question being optimized have changed.

Orient each comparison so that its original nodal contrast is positive. Each independent unpinned cell can erode that margin by up to its separately certified loss; feasible positive witnesses give a lower version of the same loss. With a fixed number of pins, selecting the largest losses minimizes the unpinned sum exactly. The witness losses prove *insufficiency* for any cheaper set by producing the opposite sign; the outer losses show *sufficiency* for the listed set. At the stated resolution these two cardinalities coincide for all five points. Units are counts of selected common cells, not equal computational cost in an evolution scheme. Tests independently exhaust small subsets to confirm the ranking logic.

The result is a common-phase, held-archive diagnostic at five imposed intervals and ordinary-floating-point enclosure precision. It does not prescribe an interacting solver grid, establish a physical transport sign or constrain cosmology.

```bash
cd kinetics
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
d=2026-11-28; p=results/early-pole-fine-charge2-scalar4; w=results/early-pole-warm
python -W error -m unittest -v test_early_pole_endpoint_warm_streaming_pin_budget
python -W error early_pole_endpoint_warm_streaming_pin_budget.py "$p-x05-state-$d.json.gz" \
 "$w-three-chart-pressure-constrained-streaming-$d.json" "$p"-gauss{8,9}warm-x05-state-$d.json.gz
```
