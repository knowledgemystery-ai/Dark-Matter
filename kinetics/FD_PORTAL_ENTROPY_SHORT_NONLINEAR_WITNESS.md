# Short nonlinear continuations confirm the original four-window signs

The [local tangent receipt](FD_PORTAL_ENTROPY_SHORT_TANGENT_WINDOWS.md) predicts three gains and one contraction on the original PCHIP plane. From the same archived x=.99 state and rebased fixed-initial-total plane, independently integrate positive and negative **finite** perturbations in the two plane directions, on the original nonlinear equations. The independent unperturbed control is integrated on the same clock; divide the central column differences by its local comoving state. The final metric is evaluated on the original local path. Use two perturbation sizes and step limits, without subtracting long .55-to-.99 histories.

| Window | Largest tangent gain minus 1 | Largest gain minus 1 estimated from finer nonlinear pair differences |
|---|---:|---:|
| .99 to hot-side reference | +3.3501643e−6 | +3.3501608e−6 |
| Across the original knot clock | +9.0257610e−8 | +9.0261465e−8 |
| Cool-side reference to last clock | −6.0454569e−7 | −6.0454881e−7 |
| Whole short window | +2.8305019e−6 | +2.8304990e−6 |

Across both directional gains and all four windows, the finer finite estimate differs from the tangent by at most 2.14e−11, and the two finite settings differ by at most 2.61e−11. This is an independent short initial-value check of this numerical tangent, not a claim about finite-size nonlinear norms of every perturbation, a new full .55-to-.99 evolution, or physical stability. The smooth cubic tangent remains a separate calculation.

```bash
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error -m unittest test_fd_portal_entropy_short_nonlinear_witness
python -W error fd_portal_entropy_short_nonlinear_witness.py "$EOS" results |
  cmp - results/fd-portal-entropy-short-nonlinear-witness-2026-11-28.json
```
