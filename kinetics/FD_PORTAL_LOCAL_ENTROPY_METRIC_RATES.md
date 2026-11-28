# Endpoint contraction does not establish monotone entropy-norm decay

The largest measured fixed-initial-total-energy endpoint gain was 0.80932762 across x=.55→1.20. Here we ask a different question: can an infinitesimal perturbation at an interior state increase the instantaneous comoving entropy norm?

For local fractional coordinates, use L=J−diag(Q′/Q), with J the fixed-basis Jacobian at the same state and Q the unperturbed comoving number, dark energy, and total energy. The squared entropy norm is eᵀGe. The tested derivative form is S=G′+LᵀG+GL; the largest generalized eigenvalue of (S,G) is the largest fractional instantaneous rate of that squared norm per unit x. G′ is differentiated along the *actual* numerical flow direction obtained from the same rate equations, not by interpolating G between sparse archive times. The EOS and positive FD bath independently recover this flow direction.

Of 17 archived combined-path samples, two give a positive full three-dimensional rate:

| x | Max full squared-norm rate | Max rate with δ total energy fixed at this same x |
|---:|---:|---:|
| .61 | +0.01235354 | −1.34177493 |
| .99 | +0.04367067 | −0.30642663 |

The other 15 full sample maxima are negative. The instantaneously fixed-total-energy plane is negative at all 17 samples, but it need not coincide with the plane transported forward from fixed initial energy at x=.55. Thus these two positive-rate directions are local counterexamples to full-space monotone contraction for the exact tested numerical metric and EOS, not counterexamples within the previously evolved initial subspace. They are not a physical instability or a claim over unsampled points. The largest combined quadrature/perturbation discrepancy over all rates is 5.88e−7, the moving-metric derivative refinement discrepancy is at most 1.91e−10, and a separate symmetric directional rate calculation agrees to 4.71e−10.

## Reproduce

From the kinetics directory, use the central normalized EOS from the previous expansion receipts (SHA-256 906d2b43bfd925c65aea3ff25ce2840cf5da7a49fcf31d8190199ae1004168ad):

```bash
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error -m unittest test_fd_portal_local_entropy_metric_rates
p=results/fd-portal
d=2026-11-28
python -W error fd_portal_local_entropy_metric_rates.py "$EOS" \
  "$p-direct-expansion-055-096-$d.json" "$p-direct-expansion-096-120-$d.json" \
  "$p-comoving-response-entropy-metric-$d.json"
```
