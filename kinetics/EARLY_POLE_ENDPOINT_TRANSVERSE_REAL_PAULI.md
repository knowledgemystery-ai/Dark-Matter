# Use real Pauli coordinates for the optional transverse spin sum

The complex factorization already avoids constructing a per-node Dirac vertex. A further algebraic reduction writes every transverse spinor bilinear as `a I + i σ·b` with real a and b. The Pauli identity

```text
(σ·u)(σ·v) = (u·v) I + i σ·(u×v)
```

reduces the scalar and vector coefficients to real dot and cross products. The Frobenius spin sum is then `2(a² + b·b)` separately for each of the two transverse polarizations. Complete the combined amplitude before applying this norm; never sum squared diagrams separately. The implementation uses the same stable finite-mass denominators and keeps the same explicit nonsingular timelike scope.

| Exact native point fixture | Complete mirrored angular mesh | Max relative coefficient difference vs complex factorization | Complex / real seconds |
|---|---|---:|---:|
| D4 at q=4.1e−6 | (4,10,18) | 3.89e−15 | 1.10 / .65 |
| D8 at q=4.35e−6 | (8,16,24) | 1.78e−15 | 15.35 / 7.40 |

Both methods use the **same** mirror and compare unchanged supports and affine forward, inverse and slope values. The timings are illustrative (roughly 1.7× and 2.1× in these two comparisons) and do not include full-band runtime. Synthetic collinear near-null and generic points have separate rounding tolerances. The original reference caller and ongoing native jobs are untouched.

```bash
cd kinetics
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
python -W error -m unittest test_finite_mass_qed_transverse_real_pauli_reference
for d in 4 8; do
  python -W error early_pole_endpoint_transverse_real_pauli_fixture.py results "$d" --timing
done
```
