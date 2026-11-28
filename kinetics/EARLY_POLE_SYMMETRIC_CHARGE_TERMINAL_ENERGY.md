# A falling scalar pressure does not imply a falling scalar energy

At the exact half-upper terminal state from the [fixed-mass pressure audit](EARLY_POLE_SYMMETRIC_CHARGE_PRESSURE.md), the retained scalar number increases and its pressure falls. Check separately whether this means that the scalar energy also fell. Use the same native quadrature as the evolution; use the nodal comparison of pressure here, not the integrated histogram reconstruction. All listed terms share archive units.

| Half-upper minus nominal at x=.05 | Energy or pressure-term difference |
|---|---:|
| Scalar energy at fixed nominal pole | +2.88501e−10 |
| Three times scalar pressure, same fixed pole | −1.32489e−11 |
| Difference of those two scalar moments | +3.01749e−10 |
| Additional scalar-energy change from the different endpoint pole | +1.00585e−10 |
| Change in the subtracted mean-field term, with its minus sign | −1.00584e−10 |
| Charged daughter energy | +5.32961e−7 |
| Charged parent energy | −5.33534e−7 |
| Sum of the charged energies | −5.72331e−10 |
| Total corrected energy | −2.83830e−10 |

The fixed-pole pressure and energy can go in opposite directions because the pressure kernel is q²/E, whereas the energy kernel also contains m²/E. At its own pole the scalar energy rises by 3.89086e−10, but this is a sector diagnostic: the mean-field correction and the net charged energy are both required to recover the closed model. In particular the pole-only scalar term and the mean-field term almost exactly cancel.

For each archived branch, rebind the 4096-step trace to its precise signed initial state and clock, then close the independently decomposed endpoint energy against the saved trace. The half-upper initial energy differs from nominal by only +3.52e−18. Its cumulative trace difference is negative at all eight archived samples and finishes at −2.83829e−10. The terminal energy contrast minus its initial contrast and accumulated trace differs by −9.39e−16, consistent with the difference of saved numerical continuity residuals. The large parent and daughter shifts do not violate the tiny total-energy contrast. This supports the internal ledger, not a separate estimate of physical thermalization, stress-energy, or spatial error.

```bash
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error -m unittest test_early_pole_symmetric_charge_terminal_energy
python -W error early_pole_symmetric_charge_terminal_energy.py results \
  results/early-pole-symmetric-charge-cutoff-audit-2026-11-28.json.gz |
  cmp - results/early-pole-symmetric-charge-terminal-energy-2026-11-28.json
```
