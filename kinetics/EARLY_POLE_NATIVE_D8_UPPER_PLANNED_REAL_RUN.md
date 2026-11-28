# Execute the fixed D8 upper plan without changing the central run

Evaluate only native full-topology IDs 157–194 from the [frozen upper raw-band plan](EARLY_POLE_NATIVE_D8_RAW_BAND_PLAN.md). The accepted [fine interior](EARLY_POLE_NATIVE_D8_UPPER_FINE_CANARY.md) and [coarse first-panel](EARLY_POLE_NATIVE_D8_UPPER_EDGE_COARSE_CHECK.md) arithmetic checks support a separate experiment using real-Pauli transverse reduction with exact azimuth reflection. Reconstruct the same state, fixed pole, topology and occupied scalar bins as the stock pilot; use the same two prescribed radial/angular meshes. The old central occupation-weighted ledgers remain stock and untouched. The stock-only upper validator also remains unchanged.

Each worker owns one independent mesh directory and alternates its planned edges toward the center. An atomic, checksummed checkpoint for each radial node binds the method and source, plan, exact native state and topology, mesh, panel, q and weight. The accepted **real-Pauli** nodes from two existing canaries seed one panel per mesh only after those canary receipts replay, the converted contributions agree with the stock references, and the full current contract is installed. No stock coefficient row enters the new ledger. The seed panels differ; their old comparison alone says nothing about same-panel angular convergence.

Reduce only when both meshes have every required node of all 38 panels. Recover raw weights by algebraically removing the held finite scalar Bose factors one exact occupied bin per panel, and check raw forward/inverse/affine as well as energy/centered-energy identities. Require at most (6	imes10^{-5}) observed two-mesh relative difference for every retained raw entry and the aggregate. If the gate fails, say underresolved; do not convert a partial or failed ledger into an accepted result. This is a numerical comparison on a fixed state and disjoint retained support, not a physical source, a total over all q, or a predicted D8 raw coverage fraction.

```bash
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error -m unittest test_early_pole_native_d8_planned_upper_real
scratch=$(mktemp -d)
python -W error early_pole_native_d8_planned_upper_real.py results "$scratch" mesh1 &
pid1=$!
python -W error early_pole_native_d8_planned_upper_real.py results "$scratch" mesh2 &
pid2=$!
wait "$pid1" && wait "$pid2" &&
  python -W error early_pole_native_d8_planned_upper_real.py results "$scratch" reduce
```

`seed1` and `seed2` independently exercise the two frozen checkpoint replays without starting new angular evaluations. No final upper-band result is included with the runner.
