# The scalar Hartree root shift largely follows the cold-count change

On the independently evolved native fine D1, D2, D4, and D8 x=.05 states, the same scalar binning, weights, bare mass, and coupling are used. Only the charge-grid resolution changes. At the old pole mass, use an exact accounting of the change in the inverse-energy gap source: cold count, normalized shape inside the same q≤1e−5 group, and all remaining modes. The cold shape term uses the regular inverse-energy deficit q²/[m E(E+m)] to avoid subtracting nearly equal terms.

Separately solve a positive scalar root after rescaling only the old cold occupation vector to the new cold count, leaving its normalized cold shape and outside modes unchanged.

| Charge step | Actual relative cold-count change | Actual relative scalar pole change | Error of count-only prediction, relative to pole **shift** |
|---|---:|---:|---:|
| D1→D2 | −1.88377% | −.63191% | 4.872e−6 |
| D2→D4 | −2.81120% | −.94599% | 4.016e−6 |
| D4→D8 | −1.37589% | −.46075% | 2.667e−6 |

The independently evaluated gap closes on all four physical snapshots. An independent bracketed root agrees with each count substitution. The derivative of the positive-gap residual is at least 2m, which gives corresponding conservative relative-to-shift bounds 7.31e−6, 6.03e−6 and 4.00e−6. For D2→D4 the normalized cold-count term is −.0281120, while the cold-shape and outside terms are −6.35e−8 and +1.75e−7. Small differences in those terms partially cancel.

This result identifies the source of the *scalar gap shift within these numerical snapshots*, not why the time-evolving charge grids change the cold count. It does not predict the full continuum D8 cold rate, which is separately measured, nor bound a physical dark-matter production rate. It also does not make the normalized cold-response profiles identical.

```bash
cd kinetics
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
python -W error -m unittest test_early_pole_endpoint_charge_cold_count_gap_ladder
python -W error early_pole_endpoint_charge_cold_count_gap_ladder.py results
```
