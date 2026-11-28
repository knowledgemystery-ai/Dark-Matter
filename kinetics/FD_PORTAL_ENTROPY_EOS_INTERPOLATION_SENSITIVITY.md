# Local entropy-norm growth is not robust to a second EOS interpolation

Keep the 17 archived visible temperatures and dark thermodynamic states fixed. Use the same normalized central entropy knots and the same low-temperature pressure anchor, but reconstruct the visible thermodynamics with a not-a-knot cubic spline in log(T) instead of PCHIP. Obtain pressure by integrating the spline-derived entropy density, then use rho=T s−P; this retains thermodynamic consistency. Checking each polynomial's extrema shows both the entropy knot function g_s and 3 g_s + d g_s/d log(T) remain positive on every retained interval. The latter is proportional to the heat capacity.

| Archived sample x | Original PCHIP maximum full squared-norm rate | Cubic-spline maximum at same held T and dark state |
|---:|---:|---:|
| .61 | +0.01235354 | −0.00705701 |
| .99 | +0.04367067 | −0.00226625 |

The two positive sampled rates both become negative; the other 15 sampled rates remain negative. All both-EOS same-time fixed-total-plane maxima are negative at these samples. The largest change in the value of the visible heat capacity among the 17 points is below 0.08%, although the time derivative of the entropy metric can change enough to reverse the local sign. The independent directional-rate consistency test is applied in the alternate interpolation as well.

This falsifies interpolation robustness of the two reported local signs, even though both remain correct for the originally specified PCHIP model. No alternate background was evolved, and this is neither a claim of monotone contraction between the samples nor a bound on physical EOS uncertainty.

```bash
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error -m unittest test_fd_portal_entropy_eos_interpolation_check
d=2026-11-28
p=results/fd-portal
python -W error fd_portal_entropy_eos_interpolation_check.py "$EOS" \
  "$p-direct-expansion-055-096-$d.json" "$p-direct-expansion-096-120-$d.json" \
  "$p-local-entropy-metric-rates-$d.json"
```

Use the same exact normalized central EOS as the original local-rate receipt.
