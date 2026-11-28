# The numerical local growth sign flips across one PCHIP temperature knot

The original x=.99 positive-rate sample has visible T=.35958348 GeV, just above the tabulated .3595 GeV entropy knot. Restart exactly at the archived .99 state and integrate the same restricted FD background only until it traverses this knot. The continuation reaches it at x=.99026086436.

The PCHIP entropy fit has a continuous value and first derivative at a knot, but its second derivative can jump. Here the exact one-sided polynomial slopes d log(c_v)/d log(T) are 3.42692203 on the hotter side and 3.36718605 on the cooler side. The consistent cubic-spline fit through the same knots gives 3.37316210 on both sides.

| Distance from crossing in x | PCHIP: hotter side | PCHIP: cooler side | Cubic: hotter side | Cubic: cooler side |
|---:|---:|---:|---:|---:|
| 8.0e−5 | +.04471511 | −.00465181 | −.00149115 | −.00131335 |
| 2.5e−5 | +.04503282 | −.00448762 | −.00125539 | −.00119982 |

Entries are largest instantaneous squared entropy-norm rates. All cubic entries use the same held states reached by the original short continuation. At each of the four states two local finite-difference settings agree to 1.11e−9 absolute or better. Derivative samples stay on one side; there is no derivative or assigned rate *at* the PCHIP knot. The sign flip is aligned with the interpolation knot and is not a physical thermal threshold. Neither fit's sign at other temperatures is certified here.

```bash
cd kinetics
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
python -W error -m unittest test_fd_portal_entropy_pchip_knot_crossing
d=2026-11-28
p=results/fd-portal
python -W error fd_portal_entropy_pchip_knot_crossing.py "$EOS" \
 "$p-direct-expansion-096-120-$d.json" "$p-entropy-eos-interpolation-check-$d.json"
```

Use the same normalized central EOS table as the earlier interpolation study.

## Independent jump audit

Only the moving visible metric coefficient appears in the jump: the normalized analytic difference of the two rate forms is −0.05317039 times the rank-one visible metric at the crossing. The smooth terms do not jump. The residual of this rank-one prediction, relative to its operator norm, falls from 0.725% to 0.226% as the offset decreases from 8e−5 to 2.5e−5. Adding this jump alone to the hot-side form yields a cold-side maximum −.00455629, within 6.87e−5 of the directly computed −.00448762 at the smaller offset. The temperature motion itself is independently reproduced from the exact short evolution.

```bash
python -W error fd_portal_entropy_knot_jump_audit.py "$EOS" \
  "results/fd-portal-direct-expansion-096-120-2026-11-28.json" \
  "results/fd-portal-entropy-pchip-knot-crossing-2026-11-28.json"
```
