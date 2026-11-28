# Re-evolving the background leaves both spline-EOS local rates negative

The previous EOS sensitivity control used a positive, thermodynamically integrated cubic spline through the same central entropy table but held the old PCHIP physical states fixed. Here independently evolve the combined FD trajectory from x=.55 to 1.20 with the spline EOS. Choose the exact same initial visible T and dark chemical/thermal coordinates; each EOS determines its own visible energy density and Hubble rate. Initial total density in the spline calculation is higher by 1.30408e−5 in relative terms.

| Sample x | Original PCHIP background and EOS | Spline EOS on held old state | Spline EOS on independently re-evolved state |
|---:|---:|---:|---:|
| .61 | +.01235354 | −.00705701 | −.00706531 |
| .99 | +.04367067 | −.00226625 | −.00227992 |

All values are largest instantaneous squared entropy-norm rates per x. The largest among all 17 re-evolved spline samples is −.00227992; no positive full-space local rate is seen at these samples. The two alternate trajectory integrations differ by at most 1.56e−10 absolute in the sampled physical coordinates. The spectrum of an audited local form is re-evaluated at the alternate physical state, not estimated by interpolating the original rate.

This strengthens the negative robustness result: changing the local thermodynamic interpolation reverses both signs, and the reversal survives re-evolution of the background. We have not shown contraction between all samples, independence from other EOS choices, or a bound on physical plasma stability. The original PCHIP results remain valid as numerical findings for that interpolation.

```bash
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
d=2026-11-28
p=results/fd-portal
python -W error fd_portal_entropy_spline_reintegration.py "$EOS" \
  "$p-direct-expansion-055-096-$d.json" "$p-direct-expansion-096-120-$d.json" \
  "$p-entropy-eos-interpolation-check-$d.json"
```

Use the same normalized central EOS table as the earlier sensitivity study.
