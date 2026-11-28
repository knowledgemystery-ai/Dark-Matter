# Cross the scalar donor and the native charged donor at x=.04

The cold-count timeline places 73.44% of the D2/D4 net terminal deficit after x=.04. The earlier frozen crossed table showed that instantaneous source changes can oppose each other. This prospective test actually continues both crossed states through the terminal interval and does not extrapolate the frozen table.

Starting from the exact native D2 and D4 x=.04 snapshots, keep each donor's own charged occupations and charge grid. Replace **only** the scalar occupation vector with the other independently evolved scalar snapshot, on their identical 67-bin scalar grid. No charged remapping is performed. Both crossed trajectories then use the same .04→.05 log-time grid and 2048-step rule as the archived native continuations. The selected fast capacity implementation is the separately checked equivalent of the bounded update and falls back to it when limiting might activate.

The setup cards verify that each scalar mass follows its scalar donor exactly at x=.04, e.g. .0006769249204 for D2 and .0006732922928 for D4. The swap deliberately changes the initial total state and energy. The energy and expansion-clock changes are recorded, not silently forced to zero. Paired number conservation, Pauli constraints, and trace closure are checked on each subsequent continuation. Do not treat the artificial mixed states or an incomplete continuation as a physical history.

The proposed endpoint comparison uses the two crossed endpoints and the two exact archived diagonal endpoints. For any terminal quantity, average each donor replacement along the two paths around the square; their sum recovers the diagonal change. Also report the finite interaction separately. Even a large effect there would be conditional on this selected .04 restart, not an attribution of the whole earlier trajectory.

```bash
cd kinetics
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
d=2026-11-28
h="results/early-pole-endpoint-charge-cold-count-history-budget-$d.json"
for c in 2 4; do
  s=$((6-c))
  python -W error early_pole_endpoint_crossed_charge_terminal_history.py results "$h" \
     "$c" "$s" "/tmp/crossed-charge-$c-scalar-$s.drive.json" --setup-only
  python -W error early_pole_endpoint_crossed_charge_terminal_history.py results "$h" \
     "$c" "$s" "/tmp/crossed-charge-$c-scalar-$s.drive.json"
done
```

Full evolved crossed endpoints are pending. Only the exact setup and protocol are recorded here.

Once both complete, the endpoint reducer takes the two crossed continuation receipts and the two exact raw final checkpoint payloads in D2-charged then D4-charged order. It independently replays the cold count, total scalar count and scalar root from each payload. It reports the native diagonal change, path-averaged donor changes and interaction separately.

```bash
python -W error -m unittest test_early_pole_endpoint_crossed_charge_terminal_receipt
python -W error early_pole_endpoint_crossed_charge_terminal_receipt.py results "$h" \
  "results/early-pole-endpoint-charge-cold-count-gap-ladder-$d.json" \
  --records C2-S4.json C4-S2.json --final-checkpoints C2-S4.drive.json C4-S2.drive.json
```
