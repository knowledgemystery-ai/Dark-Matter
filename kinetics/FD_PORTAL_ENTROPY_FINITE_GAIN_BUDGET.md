# The integrated rate terms must reproduce the finite entropy gain

Use the [validated original short tangent](FD_PORTAL_ENTROPY_SHORT_NONLINEAR_WITNESS.md). Select its single largest-gain initial direction for the whole .99-to-last window, normalize it in the initial metric, and follow that same direction throughout. At each clock evaluate the derivative of its squared entropy norm. Split the derivative into dark and visible metric motion, the visible heat-capacity log slope above three, and expansion and collision flow terms. This is an identity along a fixed direction, not an integral of independently maximized local growth scores.

Use Gauss nodes on both sides of the exact PCHIP knot with no value assigned to the knot itself. Check the density-to-thermal direction within the same smooth side; the earlier fixed-width verification stencil crossed the knot at the closest quadrature nodes and was not accepted. Compare the integrated rate with the squared endpoint gain minus one and repeat with both Jacobian and time quadrature refinements.

This attribution concerns only the numerical tangent on its held model and convention, not physical causes or a stability theorem.

```bash
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error -m unittest test_fd_portal_entropy_finite_gain_budget
python -W error fd_portal_entropy_finite_gain_budget.py "$EOS" results |
  cmp - results/fd-portal-entropy-finite-gain-budget-2026-11-28.json
```

## Completed signed budget

| Integrated term | Contribution to squared gain minus one |
|---|---:|
| Dark metric motion | +4.31691e−7 |
| Visible metric motion at reference log slope 3 | −8.49996e−5 |
| Visible heat-capacity log slope above 3 | +1.20511e−4 |
| Expansion flow | −2.40968e−5 |
| Collision flow | −6.18498e−6 |
| Sum | +5.661008e−6 |

The terms largely cancel. The derivative-sensitive visible term supplies the largest positive contribution on this trajectory; it is not the total observed amplification. The summed budget independently agrees with the direct squared endpoint gain within 3.83e−12 in both settings. This integrated convention differs from swapping EOS components at a fixed state and does not imply a physical cause assignment.
