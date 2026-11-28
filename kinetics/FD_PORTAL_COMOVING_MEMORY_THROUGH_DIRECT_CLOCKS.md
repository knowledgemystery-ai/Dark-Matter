# Initial perturbation memory through the direct expanding clock interval

Use the exact restricted combined trajectory from x=.55 to .96 and the same EOS as the archived late interval. In the original stable comoving number, dark-energy and total-energy variables, separately evolve the background and the time-ordered tangent. Compare the two columns that vary initial dark number or dark energy **at fixed initial total comoving energy** with four independent nonlinear perturbed trajectories. Total energy is not externally held fixed afterwards; pressure work remains in the evolution.

For a unit initial fractional perturbation in dark number and dark energy, the largest gain in those two final fractions at x=.96 is 1.2976283; the other singular value is about 7.30e-6 (7.46e-6 from the independent perturbations). Thus one combination is retained and slightly amplified in this normalization, while another is strongly suppressed. The ordered columns and nonlinear perturbations agree within 3.219e-7 in the initial component basis; the finest time-mesh refinements and independent kernels were compared, with no relaxed Liouville/log-determinant gate.

The independently archived .96→1.20 late operator can be composed after converting the early outputs to fractions of their background at the shared exact .96 state. The comparison is not simple multiplication of the two raw initial-normalization matrices. With the original .55 total energy fixed, the composed dark-only fractional gains at x=1.20 are 1.2387783 and about 8.49e-6, with independent local finite-perturbed-map composition giving 1.2387784 and 8.63e-6. The largest composed column discrepancy is 3.65e-7. This validates a concatenated tangent construction, not one direct nonlinear perturbation propagated across the entire concatenation.

These gains are for one dimensionless fraction norm, two selected start perturbations, and a restricted vacuum charged-lepton collision model. They are not a physical stability threshold, universal kinetic equilibrium statement or cosmological freeze-out result.

```bash
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
d=2026-11-28; p=results/fd-portal; e=eos-range.dat
python -W error fd_portal_comoving_response_through_clock_interval.py "$e" \
 "$p-direct-expansion-tracking-continuation-$d.json" "$p-direct-expansion-055-096-$d.json"
python -W error fd_portal_comoving_response_matched_composition.py "$p-comoving-response-055-096-$d.json" \
 "$p-direct-comoving-variational-propagator-$d.json"
```

Use the exact normalized public EOS input described in the tabulated FD documentation.
