# The sign of a finite PCHIP gain depends on where the window starts

The [short nonlinear check](FD_PORTAL_ENTROPY_SHORT_NONLINEAR_WITNESS.md) confirms a very small positive gain across the original knot for one fixed hot-side start. Keep the same transported plane and advance only the end of that window. Use dense short tangent solutions with two Jacobian steps and grids, and evaluate the positive endpoint metrics afresh. The gain minus one is calculated as a difference of quadratic forms, not a subtraction of nearly equal norms.

For the fixed hot-side origin x=0.9902358643642273, the largest gain minus one is about +9.026e−8 at the original cool-side endpoint. It crosses zero at approximately **x=0.990294080435**, where T≈0.3594893714 GeV. At x offsets −1e−8 and +1e−8 around that numerical crossing, the gain minus one is respectively +1.099e−10 and −1.099e−10. Both settings agree on the crossing to about 1.73e−12 in x. At the last previous clock it is about −5.152e−7.

At that same ending clock of the crossing, starting the window earlier at x=.99 instead gives a gain minus one of **+3.349e−6**. This compares the same continuous transported plane, not two unrelated bases or an arbitrary change of dynamics. It locates a finite numerical sign crossing in the stated bracket. It does not certify uniqueness, an instantaneous physical threshold, or continuous or physical stability.

```bash
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error -m unittest test_fd_portal_entropy_hot_start_finite_crossing
python -W error fd_portal_entropy_hot_start_finite_crossing.py "$EOS" results |
  cmp - results/fd-portal-entropy-hot-start-finite-crossing-2026-11-28.json
```
