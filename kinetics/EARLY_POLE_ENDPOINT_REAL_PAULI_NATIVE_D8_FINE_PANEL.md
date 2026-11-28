# Real Pauli arithmetic also reproduces the complete D8 fine shadow panel

Use the same independently completed native D8 fine panel 4, exact q support and [8,10,16,24] mesh, now with real Pauli coordinates and exact mirror pairing. All nine native affine coefficients match the original stock row within 3.22e−15 relative. The complex/mirror and real/mirror completed shadow runs took 119.78 and 61.83 wall seconds respectively; these separate illustrative runs are not a full-band runtime guarantee.

The optional real branch is explicit and limited to the tested D8 central band. It uses its own kernel source digest, disjoint from the complex factorized tag. All six wrong-method resume combinations among stock, complex and real are rejected without altering the raw line. Native current workers and their ledgers are unchanged. No D8 total is formed here.

```bash
cd kinetics
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
d=2026-11-28
p=results/early-pole
python -W error early_pole_endpoint_d8_real_pauli_fine_panel_fixture.py \
 "$p-cell-endpoint-fine-d8-x05-state-$d.json.gz" "$p-native-d8-stock-fine-panel4-$d.jsonl" \
 "$p-native-d8-real-pauli-mirror-fine-panel4-$d.jsonl"
```

Add `--compute` to recompute the real shadow. New work may explicitly select `--d8-central-cold-band --real-pauli-transverse` on a separate output; use a compatible tagged manifest and never mix methods.
