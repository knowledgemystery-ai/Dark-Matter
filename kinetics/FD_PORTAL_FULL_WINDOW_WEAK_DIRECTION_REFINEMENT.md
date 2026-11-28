# Refine the independently evolved weak FD direction

Use the same physical x=.55 seed and the published fixed-initial-total-energy weak direction. Reintegrate actual ±epsilon trajectories to x=1.20 at higher collision and bath quadrature (20, 52, 38 and 192, versus 17, 44, 32 and 168), tighter relative integration tolerance 3e−12, and max step .008. Compare both at the same epsilon and after halving epsilon under the new settings.

| Witness | Epsilon | Dark fractional column norm |
|---|---:|---:|
| Earlier direct full-window | 2e−4 | 8.12977e−6 |
| Refined direct full-window | 2e−4 | 8.30057e−6 |
| Refined direct full-window | 1e−4 | 8.56268e−6 |
| Earlier ordered tangent | zero-amplitude limit of that numerical tangent | 8.49468e−6 |

The finest tested column differs from the ordered tangent by at most 9.14e−8 absolute in any component, and its dark-only norm is about 0.80% larger. The full fractional total-energy component remains about −3.68e−5. The maximum column change at the same epsilon with new numerical controls is 2.62e−7; halving epsilon under the new controls changes the column by at most 3.75e−7. These are sensitivity tests, not an error bound or independent isolation of quadrature versus integration error. The agreement of the finest finite witness with the tangent is not proof of an asymptotic regime.

```bash
cd kinetics
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
python -W error fd_portal_full_window_weak_refinement.py "$EOS" results
```

Use the same exact normalized central EOS table as the original trajectory.
