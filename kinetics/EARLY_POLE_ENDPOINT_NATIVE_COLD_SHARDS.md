# Resume future cold panels in disjoint, verified shards

Large fixed-state cold-band quadratures are independent across native q panels. The optional shard workflow binds one exact state, pole, topology, mesh, and round-robin partition to a manifest. Each worker has a distinct native output and a local exclusive lock. A resume must find its existing panel indices inside the same shard. The join rejects missing, repeated, reordered or wrong-support rows, then restores the original native order **without reserializing numerical rows**. The ordinary one-worker, first adjacent band and topology interfaces are unchanged.

For both complete native D4 cold quadratures, split the independently archived 108 rows into four shards of 27 rows under one manifest per mesh. Replaying worker 1 on each pre-populated shard is a byte-for-byte no-op. The two independently joined full ledgers are exactly identical to their original archive bytes. Negative unit checks reject corrupted assignments, duplicate rows, and invalid topology. This establishes integrity of splitting and resumption, not an observed runtime gain or a physical convergence result.

Do not change active jobs to this mode mid-run. Use it only for new runs, or for an explicitly stopped, verified state with a separately planned disjoint allocation.

```bash
cd kinetics
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
python -W error -m unittest test_early_pole_endpoint_native_cold_shards
# STATE is the exact uncompressed state and TOP is that state's --topology-only JSON.
python -W error early_pole_endpoint_native_cold_shards.py plan "$TOP" 4 > manifest.json
for i in 0 1 2 3; do
  python -W error early_pole_endpoint_native_cold_shards.py run manifest.json "$STATE" "$i" "part-$i.jsonl" &
done
wait
python -W error early_pole_endpoint_native_cold_shards.py join manifest.json complete.jsonl part-{0,1,2,3}.jsonl
```

Run the usual complete-band coefficient and independent mesh receipts on the joined ledger before interpreting any result.

The exact archived D4 split/resume/join benchmark itself can be reproduced without repeating either integral:

```bash
python -W error early_pole_endpoint_native_cold_shard_replay_fixture.py results \
  results/early-pole-fine-d4-cold-continuum-ladder-2026-11-28.json
```
