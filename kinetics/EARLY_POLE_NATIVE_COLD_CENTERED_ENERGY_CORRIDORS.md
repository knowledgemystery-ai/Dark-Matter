# Archived centered energy widens the held D2/D4 neutral corridors

The [support-only corridor audit](EARLY_POLE_ENDPOINT_COLD_NEUTRAL_OCCUPATION_CORRIDORS.md) used complete included versus possibly intersecting numerical panels. Re-evaluate the same five subbands on both exact completed D2/D4 meshes, adding the [centered-energy partial-panel bounds](EARLY_POLE_NATIVE_COLD_SUBBAND_FROM_CENTERED_ENERGY.md). The two gross bounds remain separate; no intrapanel profile is assumed.

| q band in units of 1e−6 | Proven D2 loss and D4 gain at a common held occupation n, in units of 1e−9 |
|---|---:|
| Full [1,15] | (2.299371,2.908930) |
| [1,3] | no sufficient corridor identified |
| [3,5] | (2.387969,2.603582) |
| [5,8] | (2.586726,2.690347) |
| [8,10] | (2.747693,2.973605) |
| [10,15] | (2.920409,3.551211) |

These are open sufficient intervals, rounded inward. The central corridor is about 13.6 times wider than before; the next two bands now each have a sufficient interval. One common held occupation in (2.586726,2.603582)e−9 establishes opposite signs in the full band and in both [3,5] and [5,8]. The original diagnostic occupation 2.5862538e−9 is just below that joint interval: D2 is negative on [5,8], but D4 is still unresolved there. Full and upper intervals remain disjoint; every full-band reversal under a constant common occupation still has both upper nets negative.

This is held numerical sign information on the earlier D2/D4 charts only. It does not establish a physical occupation or a D8 sign, or modify the independent D8 comparison protocol.

```bash
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error -m unittest test_early_pole_native_cold_centered_energy_corridors
d=2026-11-28
p=results/early-pole-native-cold
python -W error early_pole_native_cold_centered_energy_corridors.py results \
  "$p-subband-from-centered-energy-$d.json" | cmp - "$p-centered-energy-corridors-$d.json"
```
