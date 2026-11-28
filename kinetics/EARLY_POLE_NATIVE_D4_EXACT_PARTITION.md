# Close the entire D4 cold number budget across five exact subbands

Combine the two separate and already admitted stock patches at (3,5,8), in units of (10^{-6}), on their original held D4 state. The original edges at 1, 10 and 15 complete the partition. The old children at 5 and 8 retain their exact geometry and coefficient bytes; their panel IDs move by one after installing the earlier cut. Require the 111 resulting native panels to be counted exactly once.

| q band (units of (10^{-6})) | Panels | Neutral occupation (units of (10^{-9})) | Held net |
|---:|---:|---:|---:|
| [1,3] | 19 | 2.23070 | (+1.52408	imes10^{-18}) |
| [3,5] | 15 | 2.30469 | (+8.68643	imes10^{-18}) |
| [5,8] | 24 | 2.43556 | (+1.40658	imes10^{-19}) |
| [8,10] | 17 | 2.64867 | (-5.17230	imes10^{-22}) |
| [10,15] | 36 | 2.92041 | (-1.14592	imes10^{-22}) |

The bandwise sum closes to the full held net within (4.63	imes10^{-33}), and total moments recover the originally unsplit ledgers within (8.44	imes10^{-15}). The central ([3,5]) band supplies **83.92% of the sum of the positive band nets**. The two upper bands together offset only **0.00610%** of that positive total. The neutral ratios are ordered as band averages, not proved to be monotonic at each momentum. A fraction of the held net is not a raw-kernel coverage fraction, a physical photon occupation, or a new D8 result.

```bash
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error -m unittest test_early_pole_native_d4_exact_partition
d=2026-11-28
python -W error early_pole_native_d4_exact_partition.py results |
  cmp - "results/early-pole-native-d4-exact-partition-$d.json"
```
