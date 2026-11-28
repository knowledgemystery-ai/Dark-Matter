# The mean-field term reverses the half-upper pressure sign

The [fixed-mass study](EARLY_POLE_SYMMETRIC_CHARGE_PRESSURE.md) found a lower scalar **kinetic** pressure at the half-upper state. That is not the full pressure of this model. The [vacuum-subtracted Hartree convention](early_pole_hartree_energy_pressure_reference.py) removes the double-counted mean-field term U from energy and adds it to pressure: ρ = Σ nE − U, P = Σ n q²/(3E) + U. The [energy replay](EARLY_POLE_SYMMETRIC_CHARGE_TERMINAL_ENERGY.md) already verifies the energy term on these exact native states.

| Half-upper minus nominal at x=.05 | Pressure contribution in archive units |
|---|---:|
| Scalar kinetic term at fixed nominal pole | −4.41631e−12 |
| Additional scalar kinetic change from using its actual pole | −6.33516e−16 |
| Hartree mean-field pressure | +1.00584e−10 |
| Corrected scalar pressure | **+9.61672e−11** |
| Corrected total modeled pressure, including both charged species | **+6.41382e−10** |

The corrected scalar pressure is up **25.57 ppm**, not down 72.5 ppm. These percentages refer to different observables and nominal denominators. The other three altered endpoints still have negative corrected scalar pressure contrasts. The half-upper **corrected total energy** remains lower by 2.83830e−10 even though its corrected total pressure rises.

The exact discrete Hartree gap removes the dynamical scalar mass from the trace: ρ−3P = Σ n m_bare²/E, using the respective bare masses for the species. Independently differentiate the closed comoving energy under a small change of log scale with occupations held fixed. Richardson-extrapolated derivatives match the same archived trace, and the differences of the paired terminal pressure and energy terms reproduce the trace contrast within 4.8e−19 in every arm.

These are the pressure and energy of the specified held numerical model. No physical plasma equation of state, trajectory stability, or spatial convergence follows.

```bash
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error -m unittest test_early_pole_symmetric_charge_hartree_pressure
d=2026-11-28
p=results/early-pole-symmetric-charge
python -W error early_pole_symmetric_charge_hartree_pressure.py results \
  "$p-cutoff-audit-$d.json.gz" "$p-terminal-energy-$d.json" |
  cmp - "$p-hartree-pressure-$d.json"
```
