# Restrict archived radial quadrature nodes without reconstructing their weights

The [centered-energy corridor](EARLY_POLE_NATIVE_COLD_CENTERED_ENERGY_CORRIDORS.md) allows a nonnegative measure anywhere inside each panel. A narrower and different question is possible on the same complete D2/D4 ledgers: what can be said about the sum of **existing Gauss–Legendre nodes** falling in each subband? The radial nodes are those of the archived quadrature (orders four and eight), not a new quadrature of the clipped panel. No individual node weight or angular sample is inferred.

For one of the positive gross coefficients, let (p) be its in-band mass fraction and let (l_i,h_i) and (l_o,h_o) be the extreme centered energies (q^2/(sqrt{m^2+q^2}+m)) of the nodes inside and outside the band. Its stored mean must lie between (p l_i + (1-p)l_o) and (p h_i + (1-p)h_o). These two inequalities yield sharp mass bounds given the means and this discrete support, independently for forward and inverse. Synthetic cases are checked against explicit linear programs. To form the net, use (S=F-B), the pointwise finite-occupation identity, so (nS-B=nF-(1+n)B). Use the same held (n=2.5862538037	imes10^{-9}), and retain the hull across meshes rather than a preferred mesh.

| q, units of (10^{-6}) | D2, node-filtered | D4, node-filtered | Continuous in-panel bounds |
|---:|:---:|:---:|:---:|
| [1,3) | unresolved | unresolved | both unresolved |
| [3,5) | negative | positive | negative / positive |
| [5,8) | negative | positive | negative / **unresolved** |
| [8,10) | negative | unresolved | negative / unresolved |
| [10,15) | negative | negative | both negative |

The new distinction is [5,8): D2's node-filtered net is in ([-1.00338	imes10^{-17},-1.71521	imes10^{-18}]); D4's is in ([7.82971	imes10^{-21},2.75117	imes10^{-19}]). The D4 neutral upper bounds on the respective discrete meshes are (2.54579	imes10^{-9}) and (2.57765	imes10^{-9}), below held (n), while the continuous-support upper bound (2.58672548	imes10^{-9}) remains slightly above it. In [8,10) the D2 masks do not split any node sets, yet D4 is still unresolved. The original continuous-support intervals remain unchanged. Certifying the continuous integral on the clipped band would require additional radial evaluation and error control; selecting existing nodes does not provide either. No inference about D8, spatial convergence, or physical photon occupation follows.

```bash
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error -m unittest test_early_pole_native_cold_gauss_node_bounds
d=2026-11-28
python -W error early_pole_native_cold_gauss_node_bounds.py results \
  "results/early-pole-native-cold-centered-energy-corridors-$d.json" |
  cmp - "results/early-pole-native-cold-gauss-node-bounds-$d.json"
```
