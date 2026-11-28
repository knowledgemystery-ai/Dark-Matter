# Most of the shifted cubic threshold lies outside the rank-one amplitude

The [crossed derivative diagnostic](FD_PORTAL_CUBIC_DERIVATIVE_CROSS.md) shows that the complete smooth-EOS threshold exceeds the PCHIP threshold by about .00202 in the near-knot layer. Separate the magnitude of the heat-capacity derivative term from everything else in that held form.

Write S(k)=B+(k−3)R. Compare both forms in common density perturbation coordinates, not with two different fractional background scales. Then the direction of R agrees numerically to 1.2e−16 in relative norm because its visible-energy covector is the same. Only its scalar magnitude differs: the cubic heat capacity is about .060% larger and the cooling magnitude about .060% smaller, making the rank amplitude about .120% smaller. Both are properties of these held model states.

| Change from PCHIP at the four held probes | Increase in critical log slope |
|---|---:|
| Replace the rank amplitude only | .0004455 to .0004483 |
| Replace the rest of the form only | .0015672 to .0015775 |
| Replace both | .0020146 to .0020277 |

When changing just the amplitude, the critical slope's **distance above 3**, not the whole slope, scales inversely with the rank amplitude. Direct crossed eigenproblems reproduce this calculation. The amplitude change explains about 22% of the shift, at the PCHIP base. It recovers a negative smooth-EOS sign on the two wider probes, but **not** on the two closer probes. The remainder of the held form is necessary to recover the latter signs. These are specific substitutions; the small interaction term means the first two changes should not be added and treated as exact independent causes.

This does not identify which individual term in the remaining form is responsible, or imply finite-window or physical stability. No rate is assigned at the PCHIP knot itself.

```bash
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error -m unittest test_fd_portal_cubic_rank_factor_attribution
python -W error fd_portal_cubic_rank_factor_attribution.py "$EOS" results |
  cmp - results/fd-portal-cubic-rank-factor-attribution-2026-11-28.json
```
