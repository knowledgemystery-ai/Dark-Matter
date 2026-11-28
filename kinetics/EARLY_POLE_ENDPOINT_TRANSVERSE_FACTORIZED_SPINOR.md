# An opt-in factorization of the transverse spin sum

The native cold continuum integrator normally constructs a four-by-four vertex for every energy and azimuth sample. For a timelike pole K=(ω,0,0,k) and either unit transverse x/y polarization ε, K·ε=0, so the slashed product anticommutes. Write the original vertex as

```text
V = [2(p·ε)/D_p + 2(b·ε)/D_b] I + [1/D_p + 1/D_b] Kslash epsslash .
```

The two complex spinor bilinears contract a constant four-by-four matrix and then combine two-by-two amplitudes before taking the absolute square. The same cancellation-avoiding denominators are retained. The optional function requires positive finite masses and μ=sqrt(ω²−k²)<2m on the outgoing branch; for future-directed on-shell b, b·K≥mμ, so the second denominator cannot cross zero on that domain. No such optimization is substituted for the separate longitudinal or audit vertex, and unmodified callers follow the original branch.

Pointwise random on-shell tests and an extreme synthetic near-null collinear case exercise reordering error separately. Actual D4/D8 paired angular probes use the exact matching state, support and quadrature, not different angular orders. The benchmark can optionally combine the algebraic factorization with exact azimuth reflection to check the composition.

```bash
cd kinetics
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
python -W error -m unittest test_finite_mass_qed_transverse_factorized_spinor_reference
python -W error early_pole_endpoint_d8_transverse_factor_fixture.py results --timing
python -W error early_pole_endpoint_d8_transverse_factor_fixture.py results --timing --mirror-factorized
```

Timing is illustrative for these point probes and independent processes on the cloud. It is not a full-band acceleration guarantee or a new D8 physical rate.

## Measured point fixtures

| Native state and q | Quadrature (cosine, inner, azimuth) | Max coefficient relative difference | Reference / optional wall seconds |
|---|---|---:|---:|
| D4, 3.5e−6 | (4,10,18) | 2.23e−16 | 5.48 / 2.81 |
| D8, 4.0e−6 | (4,10,18) | 4.45e−15 | 16.74 / 6.91 |
| D8, 4.8e−6 | (8,16,24) | 3.34e−16 | 74.85 / 40.66 |
| D8, 4.2e−6, combined with exact mirror | (4,10,18) | 2.34e−15 | 18.52 / 5.06 |

Both branches keep the same accepted angular support. These are representative point checks with the original arrays, not a change to an active production ledger. Wall timing from one run per point is illustrative only; the measured roughly 1.8–2.4× from factorization and 3.7× for the composed case should not be extrapolated to a whole band.
