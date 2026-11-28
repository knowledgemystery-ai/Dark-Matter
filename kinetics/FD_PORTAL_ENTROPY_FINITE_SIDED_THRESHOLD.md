# Switch off visible heat-capacity excess separately before or after the knot

The [uniform integrated threshold](FD_PORTAL_ENTROPY_FINITE_RANK_THRESHOLD.md) scales the excess visible metric-motion term by one multiplier. Keep the same combined FD path, same two-dimensional carried tangent and same start/end clocks, but permit two multipliers. The hot segment ends at the original PCHIP knot; the cool one starts there. No assignment of a rate at the knot is needed.

The candidate squared change in the carried Gram form is (H+alpha_h R_h+alpha_c R_c). Here (H) keeps every other integrated term, and (R_h,R_c) are the separately integrated excess terms on the two sides. Their sum must match the archived uniform calculation. Both are positive definite in the held two-dimensional frame, so increasing either coefficient cannot restore contraction. At fixed (alpha_h), contraction of every direction requires a negative definite (J=H+alpha_h R_h) and (alpha_c<1/lambda_{max}(R_c,-J)); equality gives a neutral direction. If (J) already has a positive direction, no nonnegative (alpha_c) can remove it. The same logic applies with hot/cool exchanged. This is a comparison of integrated rate allocations on an unchanged trajectory, **not** a calculation using a different EOS or a physical stability threshold.

```bash
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error -m unittest test_fd_portal_entropy_finite_sided_threshold
python -W error fd_portal_entropy_finite_sided_threshold.py "$EOS" results |
  cmp - results/fd-portal-entropy-finite-sided-threshold-2026-11-28.json
```

## Completed side-specific test

| Held excess term unchanged | Other term must fall below | Reduction of the other term at first neutral direction | Largest gain−1 when the other term is absent |
|---|---:|---:|---:|
| Hot side | Cool multiplier 0.77491152 | 22.51% of cool | (-9.66868	imes10^{-6}) |
| Cool side | Hot multiplier 0.94052995 | 5.947% of hot | (-4.13281	imes10^{-5}) |

Either side can control the sign if reduced enough in this experiment, but the same fractional change has more leverage on the hot side. The hot interval is also longer, (2.60864	imes10^{-4}) versus (8.00000	imes10^{-5}) in the expansion clock. The thresholds compare integrated scales; they do not establish how the physical EOS responds to a local change or whether the hotter interval is more effective per unit clock. Two numerical settings agree on the quoted thresholds within (5.35	imes10^{-9}). The matrix sum closes against the archived uniform calculation to numerical precision, and an independent 2×2 determinant calculation verifies the boundary.
