# The odd charge-shape contrast changes sign within the terminal history

Use the eight exactly common saved clocks of the nominal D2→D4 remap and the two accepted equal-opposite charge-shape controls. Read the stored **total** scalar yield (not a cold-cut yield at these intermediate snapshots). The raw final whole-scalar occupations independently reproduce the last entry after applying the shared normalization.

| Saved x | Half upper-minus-lower yield | Mean endpoints minus nominal yield |
|---:|---:|---:|
| .04113142 | −5.33152e−8 | −5.04480e−9 |
| .04229485 | −5.39370e−8 | −1.20594e−8 |
| .04349119 | −1.71414e−8 | −1.83892e−8 |
| .04472136 | +6.56236e−8 | −2.73545e−8 |
| .04598633 | +8.83328e−8 | −3.94772e−8 |
| .04728708 | −5.47481e−8 | −4.73427e−8 |
| .04862462 | +2.78088e−8 | −5.34001e−8 |
| .05 | +4.92226e−8 | −7.18656e−8 |

The sampled odd sequence has three adjacent sign reversals, while the even sequence is negative and moves strictly downward at all eight saved times. This is not a theorem about the continuous trajectory, higher-order response, or the total number of zero crossings. It motivates independent time refinement: the signed oscillatory part could be sensitive to numerical phase, even if a finite terminal even change remains. These observations are for one charge shape intervention, not evidence for physical production or spatial convergence.

```bash
cd kinetics
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
p=results/early-pole
d=2026-11-28
python -W error -m unittest test_early_pole_endpoint_symmetric_charge_history_budget
python -W error early_pole_endpoint_symmetric_charge_history_budget.py \
 "$p-charge-symmetric-null-terminal-$d.json" "$p-charge-symmetric-null-lower-x05-receipt-$d.json" \
 "$p-two-moment-charge-from2-x05-receipt-$d.json" "$p-charge-symmetric-null-upper-x05-receipt-$d.json"
```
