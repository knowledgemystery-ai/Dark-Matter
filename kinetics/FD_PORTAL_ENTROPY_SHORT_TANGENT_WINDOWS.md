# A local tangent resolves finite windows that long-history subtraction cannot

The [complete direct finite-difference control](FD_PORTAL_ENTROPY_FINITE_WINDOWS.md) did not resolve the sign of tiny gains across these four short intervals. Instead start from the [archived transported fixed-initial-total plane](FD_PORTAL_TRANSPORTED_ENTROPY_PLANE.md) at x=.99, on the original PCHIP state. Rebase its ill-conditioned fractional columns by QR without changing the plane. Propagate both columns with the locally evaluated fractional Jacobian, without starting again at x=.55 or projecting at either endpoint.

Use the original thermal-knot path. Split the integration at the knot but do not use a metric derivative or a finite window endpoint there. Use four and six Jacobian nodes per short segment, with two independent centered Jacobian steps, then compare the two resulting gains. If X and Y are the two endpoint entropy forms in the same transported coordinates, solve the generalized eigenproblem (Y−X,X) and recover the gain minus one as λ/(sqrt(1+λ)+1). This avoids subtracting two numbers nearly equal to one.

This is a numerical tangent of the original replayed state initialized by the archived carried plane, not a new nonlinear history, continuous monotonicity or physical stability theorem. The smooth cubic is a separate comparison.

```bash
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error -m unittest test_fd_portal_entropy_short_tangent_windows
python -W error fd_portal_entropy_short_tangent_windows.py "$EOS" results |
  cmp - results/fd-portal-entropy-short-tangent-windows-2026-11-28.json
```

## Four measured clock windows

| Window | Original PCHIP gain minus 1 | Separately continued smooth cubic gain minus 1 |
|---|---:|---:|
| .99 to hot-side reference | +3.35016e−6 | −2.23945e−6 |
| Across original knot clock | +9.02576e−8 | −4.59101e−7 |
| Cool-side reference to last clock | −6.04546e−7 | −5.06380e−7 |
| Whole short window | +2.83050e−6 | −3.20493e−6 |

The original has growth in three measured windows, including a small gain across the original knot, and contraction in the cool-side window. Its two tested settings differ by at most 1.94e−14. The [smooth cubic comparison](FD_PORTAL_ENTROPY_CUBIC_SHORT_TANGENT.md) instead initializes from its own archived plane and contracts in all four windows. Even with these finite differences in the held numerical model, no result here establishes continuous or physical stability.

**Independent check:** [Local nonlinear continuations](FD_PORTAL_ENTROPY_SHORT_NONLINEAR_WITNESS.md) reproduce all four original PCHIP signs from the same archived .99 plane.
