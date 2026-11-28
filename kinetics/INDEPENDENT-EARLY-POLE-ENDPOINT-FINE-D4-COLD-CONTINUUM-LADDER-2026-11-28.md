# The integrated cold-band contrast does not converge under D1→D2→D4 charge refinement

Integrate the common retained q interval [1e−6, 1.5e−5] at x=.05 and k=.0015. Keep the scalar grid unchanged; evolve each fine charge chart independently and use its own scalar mass, charge occupations, and matched transverse pole. On D4 the complete 108 native panels on both quadratures reproduce the same support and pole. The largest relative difference among the nine integrated coefficients is 6.93e−14 between the two D4 quadratures. That is a fixed-state integration check, not a spatial error bound.

| Charge subcells per archived cell | Complete panels | Gross forward number / original fine archive | Conditional number-neutral photon occupation | Relative split of energy- and number-neutral roots |
|---:|---:|---:|---:|---:|
| 1 | 27 | 6.45841084 | 3.34922719e−9 | +4.47663e−8 |
| 2 | 49 | 1.30271123 | 2.90893020e−9 | +5.34301e−8 |
| 4 | 108 | .07649859 | 2.29937066e−9 | +7.67485e−8 |

The D4 gross forward number is **5.872% of D2**, while the gross inverse number is **4.642% of D2**. The coarse cold-band gain relative to the old fixed archive thus changes from above unity on D1 and D2 to far below unity on D4. The selected two-point aperture result was not just a single-q outlier: the large shift survives integration on the same retained band.

For the held diagnostic photon occupation 2.58625380e−9, selected as the geometric midpoint of the D2 and D4 number-neutral roots, the net number source divided by gross inverse is −11.09% on D2 and +12.48% on D4. This is a conditional sign test at a fixed, constructed occupation, not a determination of the live photon population. The forward median enclosures for D2 and D4 overlap, so the ordering of the medians is not resolved. Above q=1e−5 the D4 fraction of gross forward number is 9.505e−6, about 1.95 times the D2 fraction; nevertheless its absolute forward tail weight is only 11.45% of D2. There is no spatial convergence, extrapolation outside the retained interval, or physical production claim.

## Reproduce

```bash
cd kinetics
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
python -W error -m unittest test_early_pole_endpoint_cold_band_panel_contrast
python -W error early_pole_endpoint_fine_d4_cold_continuum_ladder.py results \
  results/fine-d4-endpoint-continuum-cold-{1,2}.jsonl.gz
```

To recompute the D4 panels themselves, decompress the existing independently evolved x=.05 state and run both explicit integration meshes. Use new output paths unless continuing an exact matching run.

```bash
d=2026-11-28
gzip -dc "results/early-pole-cell-endpoint-fine-d4-x05-state-$d.json.gz" > /tmp/fine-d4-x05.json
for m in 1 2; do
  python -W error early_pole_endpoint_continuum_cold_band_pilot.py \
    /tmp/fine-d4-x05.json "$m" "/tmp/fine-d4-cold-$m.jsonl"
done
```
