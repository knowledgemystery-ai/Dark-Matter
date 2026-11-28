# Coordinate dependence of the same FD memory response

For the exact direct combined path and the two perturbations at fixed initial total comoving energy, pull the earlier fractional number/energy response back to variations (eta_dark, log T_dark) at each endpoint. The chart transformation is local, using the analytic number and energy susceptibilities at the exact archived states. Fourth-order independent finite differences check all three susceptibility matrices (maximum absolute defect less than 3.8e-11).

If J(x) maps intensive dark perturbations to fractional number and dark energy at x, the intensive endpoint matrix is J(end)^(-1) F J(start), where F is the dark-only fractional response. For the same x=.55 start the largest gains in these two coordinate norms are:

| Endpoint | Fractional number/energy norm | eta/log-temperature norm |
|---:|---:|---:|
| .96 | 1.2976283 | 2.6027643 |
| 1.20 | 1.2387783 | 2.5321293 |

The small intensive gains are 3.58e-6 and 4.09e-6; the independently perturbed interval maps give 3.66e-6 and 4.16e-6. Thus a strong difference between two directions survives this particular change of variables, but the numerical gain is not an invariant amplification factor. This does not convert the benchmark into a physical instability or freeze-out criterion.

```bash
cd kinetics
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
d=2026-11-28; p=results/fd-portal; e=eos-range.dat
python -W error fd_portal_composed_memory_coordinate_check.py "$e" \
 "$p-direct-expansion-055-096-$d.json" "$p-direct-expansion-096-120-$d.json" \
 "$p-comoving-response-055-096-$d.json" "$p-direct-comoving-variational-propagator-$d.json" \
 "$p-comoving-matched-response-055-120-$d.json"
```
