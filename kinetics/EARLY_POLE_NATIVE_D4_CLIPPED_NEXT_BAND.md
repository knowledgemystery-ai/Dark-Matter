# Exact-cut follow-up on the held D4 cold band [5,8] × 10⁻⁶

The earlier [node-filtered proof](EARLY_POLE_NATIVE_COLD_GAUSS_NODE_BOUNDS.md) establishes a positive D4 sum when only the **existing nodes** in this band are counted. The [continuous-support moment bounds](EARLY_POLE_NATIVE_COLD_CENTERED_ENERGY_CORRIDORS.md) permit zero at the same held photon occupation. Those statements use different information and do not conflict. To check an actual quadrature of the clipped band, install its endpoints on the native D4 graph and re-integrate the affected children. The exact admitted D4 snapshot, scalar grid, stock angular kernel, diagnostic photon occupation and all uncut rows are unchanged.

The archived parent panels 32 and 55 have bounds ([4.778274331905574,5.083582387749736]	imes10^{-6}) and ([7.863047969834538,8.046092007084980]	imes10^{-6}). Split each at the requested cut. The resulting topology contains 110 panels; the only recalculated IDs are 32, 33, 56 and 57. Both children of each parent are evaluated, including the portions outside the target span. Compare their sum to the original parent for all nine stored entries, with a relative tolerance of (6	imes10^{-5}) applied entry by entry so that small inverse entries cannot hide behind the gross forward scale. Require the same limit when comparing each retained entry across the two full angular/radial mesh choices.

For the exact-cut band, report the two observed numerical estimates and a conservative coefficient box across their observed spread: take the lower slope and upper inverse separately per panel, and vice versa for the upper net. If the box crosses zero, the sign remains unresolved. This is an observed two-mesh comparison on a fixed numerical state, **not a bound on true continuum integration error** or evidence about D8 or physical photon occupation.

```bash
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error -m unittest test_early_pole_native_d4_clipped_next_band
scratch=$(mktemp -d)
python -W error early_pole_native_d4_clipped_next_band.py results "$scratch" plan
python -W error early_pole_native_d4_clipped_next_band.py results "$scratch" mesh1
python -W error early_pole_native_d4_clipped_next_band.py results "$scratch" mesh2
python -W error early_pole_native_d4_clipped_next_band.py results "$scratch" finish
```

## Why another quadrature is necessary

On each of the two original cut panels, put smooth, normalized bumps in four short windows adjacent to the two ends and the cut. A window of width (w=10^{-10}) uses (30t^2(1-t)^2/w), (tin[0,1]); this joins zero with a continuous first derivative. For each positive slope and inverse measure, take a convex mixture of two windows that bracket its saved mean centered energy. The exact weight is the linear interpolation in the integral mean of (q^2/(sqrt{m^2+q^2}+m)). Reconstruct the forward density as the sum of slope and inverse. Total mass, centered energy, and ordinary energy then equal the archived values to closure tolerance. For a wholly contained or disjoint panel, its strictly internal archived mean can likewise be realized by two smooth interior bumps, and selecting the band contributes its whole mass or zero.

With one mixture choice the exact-band net is (-3.85266	imes10^{-22}); reversing the choices gives (+2.84120	imes10^{-19}). Both signs occur on both archived meshes with the same moments and held occupation. These are feasibility witnesses, not candidate microscopic collision kernels, and do not determine the new quadrature's sign. In particular one cannot prove the continuous-support sign from the saved moments alone simply by excluding point masses.

```bash
python -W error -m unittest test_early_pole_native_d4_cut_moment_witness
python -W error early_pole_native_d4_cut_moment_witness.py results "$(mktemp -d)" |
  cmp - results/early-pole-native-d4-cut-moment-witness-2026-11-28.json
```

As a separate sensitivity diagnostic, holding the upper endpoint fixed gives the zero of the original conservative D4 lower bound at a *different* lower cutoff, (q=5.00069096414	imes10^{-6}). On this local branch the lower panel can place zero slope above the cut, while its maximal inverse tail falls inversely with the difference in centered energies between the cut and its lower edge. Evaluating the general moment enclosure confirms the analytic root on both meshes. The slight aperture change is not a substitute for the specified cut.

## Completed quadrature on the original aperture

Both split meshes retain 24 panels inside exactly (q=[5,8]	imes10^{-6}). The nine entries reclose their unsplit parents within (3.56	imes10^{-14}) componentwise and agree across the two meshes within (2.72	imes10^{-13}) on each retained panel. The observed per-panel two-mesh coefficient box gives (nS-B=[1.4065781567228356,1.4065781567287908]	imes10^{-19}), strictly positive on the held D4 state at the original test occupation. Earlier moment-only bounds already put D2 below zero; the independent D4 re-integration now has the positive numerical sign at the same aperture. This does not promote an observed two-mesh box to a rigorous integration error bound, a D8 comparison, or a claim about the physical occupation.

```bash
d=2026-11-28; scratch=$(mktemp -d)
cp "results/d4-clipped-next-band-patches-$d/"*.jsonl "$scratch/"
python -W error early_pole_native_d4_clipped_next_band.py results "$scratch" finish |
  cmp - "results/early-pole-native-d4-clipped-next-band-$d.json"
```
