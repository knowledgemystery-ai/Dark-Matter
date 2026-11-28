# Visible entropy-metric motion raises the cubic threshold; expansion partially offsets it

The [rank-amplitude split](FD_PORTAL_CUBIC_RANK_FACTOR_ATTRIBUTION.md) separated the heat-capacity derivative coefficient from the rest of the local rate form, B=S(k=3). Decompose B into four additive terms, retaining common density perturbation coordinates and the original PCHIP rank amplitude. The dark and visible terms come from the derivative of the entropy metric; the expansion and collision terms are the corresponding fractional flow Jacobians contracted with the local metric. The visible term excludes the change from k=3 already singled out in the rank term.

At each of the four held states near the knot, evaluate all 16 subsets that swap these terms from PCHIP to cubic. Every mixed reference form stays strictly negative at k=3, so the same rank-one threshold formula applies. Average each term's increment over all replacement orders. The allocations close to the independently reproduced threshold difference. Individual substitutions at the PCHIP base give the same dominant pattern.

| Term replaced | Allocated shift in critical log slope |
|---|---:|
| Visible metric motion, excluding k−3 | +.0028108 to +.0028297 |
| Expansion flow | −.0012523 to −.0012436 |
| Collision flow | +6.96e−8 to +6.99e−8 |
| Dark metric motion | +4.48e−10 to +7.84e−10 |

The net shift of .00157 outside the rank amplitude is mostly the **time variation of the visible entropy metric**, partially cancelled by the expansion part of the held flow. The direct collision-flow swap barely moves this threshold in these particular comparisons. This is not a bound on collision strengths, all perturbation modes, or the contribution of collisions to the past trajectories.

All comparisons use the same held states and an algebraic replacement convention. They are not new EOS histories, a determination of physical stability, or a substitute for finite-time perturbations. No rate is assigned at the interpolation knot.

```bash
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error -m unittest test_fd_portal_cubic_remainder_attribution
python -W error fd_portal_cubic_remainder_attribution.py "$EOS" results |
  cmp - results/fd-portal-cubic-remainder-attribution-2026-11-28.json
```
