# The cold-count scalar root substitution also preserves the tested raw aperture

At the exact two archived native D4 and D8 x=.05 endpoints, retain the late charged occupations and matched transverse photon pole in both evaluations. Replace only the scalar mass by the earlier cold-count-only Hartree prediction from the adjacent lower chart. Compare the raw forward, inverse and affine-slope traces at fixed k=.0015 and at q=3.5e−6 and 8.2e−6.

| Chart | q | Relative forward difference | Relative inverse difference |
|---:|---:|---:|---:|
| D4 | 3.5e−6 | −1.27849e−5 | −1.38553e−5 |
| D4 | 8.2e−6 | −1.13855e−5 | −1.23034e−5 |
| D8 | 3.5e−6 | −4.60759e−6 | −5.07301e−6 |
| D8 | 8.2e−6 | −4.03940e−6 | −4.38806e−6 |

The two numerical quadratures for each arm differ by less than 5.62e−13 relative. The count-based mass is therefore sufficient for these **four tested raw point apertures** when the other state variables are held fixed. This does not remove the charged-history differences, certify an integral over q, or bound the source with a physical scalar Bose population. The complete D8 central-band run is a separate calculation.

```bash
cd kinetics
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
d=2026-11-28
python -W error early_pole_endpoint_charge_cold_count_aperture.py results \
  "results/early-pole-endpoint-charge-cold-count-gap-ladder-$d.json" \
  "results/early-pole-fine-d4-cold-continuum-ladder-$d.json"
```
