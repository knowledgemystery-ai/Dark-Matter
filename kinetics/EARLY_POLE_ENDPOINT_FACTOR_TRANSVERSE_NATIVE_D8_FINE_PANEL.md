# A complete native D8 fine panel also replays at unchanged quadrature

As a separate shadow test, choose already completed stock panel **4** of the high native D8 mesh [8,10,16,24]. Its exact q window is [3.193617887418766e−6, 3.287067390627979e−6]. Evaluate the same complete panel with the existing opt-in complex factorization and exact azimuth pairing. This is not a different q, radial, cosine or azimuth integration order.

The nine full affine coefficients match to 1.12e−15 relative. The original line and independently computed alternative are stored exactly with separate source/kernel tags. The complete stock run continues without replacement and its incomplete wider fine ledger is **not** integrated or interpreted here. One completed fine panel validates this code path, not the entire band or a speed guarantee.

```bash
cd kinetics
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
d=2026-11-28
p=results/early-pole
python -W error early_pole_endpoint_d8_factorized_fine_panel_fixture.py \
 "$p-cell-endpoint-fine-d8-x05-state-$d.json.gz" "$p-native-d8-stock-fine-panel4-$d.jsonl" \
 "$p-native-d8-factorized-mirror-fine-panel4-$d.jsonl"
```

Add `--compute` to reproduce the optional panel from scratch, rather than checking the stored raw lines.
