# One complete native D8 panel replays with optional factorization and azimuth reflection

Lock already-completed stock coarse panel **10** of the native D8 central band, on [3.547360292955679e−6, 3.551692092529729e−6] and mesh [4,6,10,18]. Independently recompute it into a separate shadow file, with transverse factorization and the analytically paired azimuth samples. This changes neither the native q topology nor the quadrature order. The nine accumulated affine number, energy, and centered-energy coefficients agree to at most 1.23e−15 relative against the nominal calculation. The shadow bytes are exactly reproducible from the same inputs.

The serial and sharded interfaces recognize the factorized kernel only when explicitly enabled on the native D8 central band. Output rows bind the factorized kernel source digest and the mirror flag. A stock worker cannot reuse a tagged row, and a factorized worker cannot reuse a stock row; a correctly tagged, complete shard resumes unchanged. The original ongoing D8 job is not interrupted and no stock results are replaced. This one-panel check does not provide the D8-band integral or a full-run speed guarantee.

```bash
cd kinetics
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
p=results/early-pole
d=2026-11-28
python -W error early_pole_endpoint_d8_factorized_panel_fixture.py \
 "$p-cell-endpoint-fine-d8-x05-state-$d.json.gz" "$p-native-d8-stock-panel10-$d.jsonl" \
 "$p-native-d8-factorized-mirror-panel10-$d.jsonl" --compute
```

For normal future runs use `--d8-central-cold-band --factorized-transverse` on a separate output; use `--mirror-azimuth` only when intentionally selecting the corresponding tagged method in the exact manifest. Never join physical ledgers from different kernels or mirror settings.
