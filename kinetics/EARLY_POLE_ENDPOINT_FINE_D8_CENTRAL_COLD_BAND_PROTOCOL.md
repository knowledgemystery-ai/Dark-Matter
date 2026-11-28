# Prespecified D8 cold-band follow-up

The completed D2 and D4 results enclose more than 80% of each chart's retained gross forward number contribution inside q=[3e−6,5e−6]. On that measured basis, choose this common central band before evaluating integrated D8 coefficients. It is a targeted next-resolution check, not the full cold band [1e−6,1.5e−5].

The exact independently evolved native fine D8 x=.05 state (2048 steps, unchanged scalar grid) yields **34** panels on the selected band, compared with 195 on the full D8 cold span. A 49/89 internal topology comparison and both published integration-mesh topology cards agree on all shared panel bounds, the matched pole, and the D8 endpoint fingerprint. The special-band option accepts only this native D8 chart geometry and versioned step count, stores a canonical state digest on every row, and keeps the existing main and D2 adjacent-band interfaces. The total forward/inverse coefficients and any conditional sign remain **pending** until both mesh ledgers pass the complete-panel validator. Its D4 reference is the previous enclosure over the same q interval, without inventing D4 intrapanel shapes.

```bash
cd kinetics
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
python -W error -m unittest test_early_pole_endpoint_fine_d8_central_cold_band
d=2026-11-28
p=results
gzip -dc "$p/early-pole-cell-endpoint-fine-d8-x05-state-$d.json.gz" > /tmp/fine-d8-x05.json
for m in 1 2; do
  python -W error early_pole_endpoint_continuum_cold_band_pilot.py /tmp/fine-d8-x05.json "$m" \
    --d8-central-cold-band --topology-only
  python -W error early_pole_endpoint_continuum_cold_band_pilot.py /tmp/fine-d8-x05.json "$m" \
    "/tmp/fine-d8-central-cold-$m.jsonl" --d8-central-cold-band
done
python -W error early_pole_endpoint_fine_d8_central_cold_receipt.py /tmp/fine-d8-x05.json \
  "$p/early-pole-fine-d4-cold-subband-enclosure-$d.json" /tmp/fine-d8-central-cold-{1,2}.jsonl
```

Do not treat topology, a partial ledger, or the relative stability of the original fixed-state quadratures as a D8 convergence result. The held occupation in the eventual sign test remains a construction from D2/D4, not a live photon population.
