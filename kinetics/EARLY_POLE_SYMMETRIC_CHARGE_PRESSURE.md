# More scalar particles can coexist with lower terminal fixed-mass kinetic pressure

Use the exact five matched 4096-step states already replayed by the [terminal cutoff audit](EARLY_POLE_SYMMETRIC_CHARGE_CUTOFF_AUDIT.md), at x=.05 and fixed nominal measurement mass 0.0008674814895883813. This is a pressure moment of the **saved scalar component**. At the common x, the overall dimensionful conversion and the factor of one third cancel in relative comparisons.

| Arm | Whole retained count vs nominal | Fixed-mass pressure vs nominal (ppm) | Pressure interval with independent shapes and fixed node mean energies (ppm) |
|---|---:|---:|---:|
| Full lower | −332.249 ppm | −432.025 | [−432.151, −431.913] |
| Full upper | −62.157 ppm | −337.530 | [−337.669, −337.431] |
| Half lower | −118.862 ppm | −119.780 | [−119.899, −119.661] |
| Half upper | +20.385 ppm | −72.515 | [−72.641, −72.403] |

The central pressure column assumes constant occupation per native cell and independently re-integrates the cell average at 16 and 32 Gauss nodes. The four relative contrasts agree within 5.5e−20 absolute. In the half-upper run the pressure per particle is down 92.899 ppm. Its cold region q≤1e−5 contributes a small **positive** 0.0268 ppm to the total pressure difference. The largest negative one-cell contribution comes from q=[2.74240,3.96368], at −24.064 ppm. Replacing the held mass only by that branch's saved mass at fixed nodal occupations gives −0.0104 ppm: the measured contrast is not driven by that mass change.

The half-upper sign survives either of two reconstruction conditions, not unrestricted independent shape:
- If the arms share a conditional shape in each cell, multiply the signed cell count differences by the endpoint values of the increasing kernel q²/E. The half-upper interval is [−91.459,−50.468] ppm relative to nominal histogram pressure.
- If each arm instead redistributes its cell mass independently, hold its mean energy in that cell at the native node energy. The pressure kernel as a function of energy, E−m²/E, is concave, so its mean lies between the edge chord at the fixed mean energy and the function value at the mean. Subtract the two arm bounds crosswise. The last column is then relative to nominal **nodal** pressure.
- If only the cell counts are fixed and both arms can have arbitrary independent shapes, the broader bounds admit both pressure signs. The number result alone is not a pressure certificate.

These are conditional moments of completed numerical trajectories, not a physical pressure estimate, whole-system stress-energy test, or spatial convergence result.
The [self-consistent Hartree pressure follow-up](EARLY_POLE_SYMMETRIC_CHARGE_HARTREE_PRESSURE.md) adds the mean-field term and reverses the half-upper scalar pressure sign. Do not use this fixed-mass diagnostic as the pressure of the closed model.

```bash
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error -m unittest test_early_pole_symmetric_charge_pressure
python -W error early_pole_symmetric_charge_pressure.py results \
  results/early-pole-symmetric-charge-cutoff-audit-2026-11-28.json.gz |
  cmp - results/early-pole-symmetric-charge-pressure-2026-11-28.json
```
