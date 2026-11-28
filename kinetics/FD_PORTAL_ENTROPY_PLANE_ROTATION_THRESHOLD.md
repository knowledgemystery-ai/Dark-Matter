# The plane carried from the initial state is not the instantaneous energy plane

The [transported-plane experiment](FD_PORTAL_TRANSPORTED_ENTROPY_PLANE.md) found a positive local entropy-norm rate at x=.99 when fixing **initial** total energy, not when fixing the **current** total energy. An invariant angle in the positive entropy metric quantifies how different these two restrictions are. This test does not change the archived flow or perturb either plane.

| Held x | Smallest angle needed for any nonnegative-rate direction | Actual carried vs current plane, same entropy metric | Largest actual carried-plane rate |
|---:|---:|---:|---:|
| .61 | [84.49111°, 84.49181°] | 2.78753° | −1.341773 |
| .99 | [68.48639°, 68.48711°] | 82.72803° | +.027685 |

Whiten the three fractional perturbations by the square root of the entropy metric. The symmetric rate then has exactly one positive eigenvalue at each held point. The current plane is strictly negative. Each connected half of the nonnegative-rate cone therefore stays on one side of that plane. Its closest unit direction must be on the cone boundary: the only positive interior stationary value of the signed unit normal dot product would be its maximum.

In the eigenbasis, the zero-rate boundary is the normalized ellipse (a cos φ, b sin φ, 1), with a and b given by square roots of the positive eigenvalue divided by the respective absolute negative eigenvalues. The dot product with a unit plane normal has Lipschitz constant at most max(a,b), so the nearest sample within π/N bounds the possible missed minimum. These angle intervals refer to the held matrices; they do not bound their continuum discretization error.

At .61 the carried plane is too close to the strictly negative plane to contain any nonnegative direction. At .99 it has moved farther than the minimum required, and its measured orientation actually intersects the growth cone. Exceeding the minimum alone would not guarantee growth for an arbitrary orientation. The two angles compare local constraints, not a physical rotation associated with observable instability or the finite-time propagator.

```bash
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error -m unittest test_fd_portal_entropy_plane_rotation_threshold
python -W error fd_portal_entropy_plane_rotation_threshold.py "$EOS" results |
  cmp - results/fd-portal-entropy-plane-rotation-threshold-2026-11-28.json
```
