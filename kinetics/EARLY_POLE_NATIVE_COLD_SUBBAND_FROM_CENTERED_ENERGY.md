# The centered-energy panels resolve opposite D2 and D4 net signs in the central band

The [positive-panel enclosure](EARLY_POLE_ENDPOINT_FINE_D4_COLD_SUBBAND_ENCLOSURE.md) used only whole included and possibly intersecting panels. It found the central D4 net positive but left the central D2 sign unresolved. Add the stored centered-energy moments already used to [bound the mean momentum](EARLY_POLE_NATIVE_COLD_MEAN_FROM_CENTERED_ENERGY.md).

Within a crossing panel, let [u_l,u_r] be the exact range of u=E−m, u_c its value at the subband cut, and μ the normalized centered-energy moment. If probability p lies below the cut, limiting mean values are at least p u_l+(1−p)u_c and at most p u_c+(1−p)u_r. Invert these inequalities and clip to [0,1]. This gives more information about partial panel mass without an intrapanel density fit. Enclose forward, inverse and stimulated slope separately using their own positive moments.

| Native chart | Weighted forward enclosure, q=[3e−6,5e−6] | Conditional net enclosure at the original held occupation |
|---|---:|---:|
| D2 | [4.98265,5.58352]e−7 | [−3.36904e−16,−9.67559e−18] |
| D4 | [2.99876,3.10760]e−8 | [+5.94610e−18,+1.12124e−17] |

Both complete accepted meshes give D2 negative and D4 positive at the **same previously held** n=2.5862538e−9. The earlier support-only D2 upper bound was positive; its sign was not assignable without additional moments. These are tighter numerical enclosures on the selected held D2/D4 charts, not a D8 result, a change to the original D8 admission gate, or a physical photon sign prediction.

```bash
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error -m unittest test_early_pole_native_cold_subband_from_centered_energy
d=2026-11-28
p=results/early-pole-native-cold
python -W error early_pole_native_cold_subband_from_centered_energy.py results \
  "$p-mean-from-centered-energy-$d.json" | cmp - "$p-subband-from-centered-energy-$d.json"
```
