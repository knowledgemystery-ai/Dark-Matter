# Test exact real-Pauli arithmetic on one disjoint upper D8 panel

The [D8 upper raw-band plan](EARLY_POLE_NATIVE_D8_RAW_BAND_PLAN.md) is separate from the still-running, occupation-weighted [D8 central band](EARLY_POLE_ENDPOINT_FINE_D8_CENTRAL_COLD_BAND_PROTOCOL.md). Before a wider upper-band calculation, take the single full-topology panel containing the independently published upper test momentum (q=1.3344923953650561	imes10^{-5}). This selects index 183 with edges ([1.3269369822998557,1.3592223282582193]	imes10^{-5}) on the unchanged native (x=.05), 2048-step D8 state.

On the finest archived angular/radial setting ((8,10,16,24)), evaluate the original stock implementation and separately evaluate the real-Pauli transverse reduction with exact azimuth reflection. Preserve the stock row and all eight real-Pauli nodes, including their state digest, method digest, bounds and weights. Require relative agreement below (2	imes10^{-12}) for **all nine** integrated entries (gross forward, inverse and affine slope; each with number, energy and centered-energy weights). Validate the exact full topology and pole independently, and reject mixed node provenance. The original central ledger is neither read nor altered.

```bash
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error -m unittest test_early_pole_native_d8_upper_fine_canary
d=2026-11-28; scratch=$(mktemp -d)
# Replay held exact checkpoints without repeating the expensive angular work:
cp "results/d8-upper-fine-canary-checkpoints-$d/"* "$scratch/"
python -W error early_pole_native_d8_upper_fine_canary.py results "$scratch" |
  cmp - "results/early-pole-native-d8-upper-fine-canary-$d.json"
```

Starting with an empty scratch directory reruns the ambient stock and real-Pauli quadratures. This is an arithmetic comparison on one physical-kernel panel, not a second angular-resolution rung, a total or raw coverage on the upper band, spatial convergence, or a physical photon history.

## Completed panel comparison

The largest componentwise relative discrepancy was (1.11	imes10^{-15}). The held weighted stock gross forward and inverse number entries were (1.07931	imes10^{-14}) and (2.88770	imes10^{-23}), respectively, for this one panel. Independently reconstructing the nine integrated entries from the explicit node ledger reproduces both implementations. A changed node momentum is rejected on checkpoint replay. No wider upper-band computation or modification of the weighted central run is part of this artifact.
