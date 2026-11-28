# The separately evolved cubic plane contracts over all four short windows

The [direct finite-history control](FD_PORTAL_ENTROPY_FINITE_WINDOWS.md) was too imprecise to resolve short-window ratios. Start instead at x=.96 with the archived plane in the [smooth cubic endpoint tangent](FD_PORTAL_SPLINE_ENTROPY_ENDPOINT_GAIN.md). Rebase the two columns and re-evolve only the short remaining interval in the **actual smooth cubic EOS**, matching its own saved x=.99 state. Do not borrow the PCHIP plane or reinterpret the cubic as having a knot at the original crossing time.

| Same held clock window as prior comparison | Largest cubic gain minus 1 |
|---|---:|
| .99 to hot-side reference | −2.23945e−6 |
| Across the original PCHIP knot clock | −4.59101e−7 |
| Cool-side reference to last clock | −5.06380e−7 |
| Whole short window | −3.20493e−6 |

Two distinct background ODE settings, Jacobian steps and sample meshes agree within 6.8e−15 on the reported tangent gains. The alternate held plane contracts in all four measured windows. This is not a new full nonlinear trajectory or a statement of physical stability or monotonicity on any different span.

```bash
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error -m unittest test_fd_portal_entropy_cubic_short_tangent
python -W error fd_portal_entropy_cubic_short_tangent.py "$EOS" results |
  cmp - results/fd-portal-entropy-cubic-short-tangent-2026-11-28.json
```

**Independent check:** [Short cubic nonlinear continuations](FD_PORTAL_ENTROPY_CUBIC_SHORT_NONLINEAR.md) reproduce all four signs on its own alternate .99 state.
