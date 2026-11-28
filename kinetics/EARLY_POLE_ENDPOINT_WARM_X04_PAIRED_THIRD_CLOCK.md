# Actual matched 4096-step contrast for G9 and mixed300

Repeat the exact accepted mixed300 restart from x=.02 through x=.04 at 4096 uniform log-time steps and the original fixed one-cell upgrade. This is independently evolved, not obtained by applying an interpolation to the G9 result. The older mixed300 and G9 results at 1024 and 2048 use their respective unchanged restart states. All reported pressures use the same native measurement mass.

For 2048 divided by 4096 minus one, the number shifts are −3.17156958e-5 (G9) and −3.17156952e-5 (mixed300). By contrast, the matched mixed300/G9 count ratio changes by only −5.84643e-13 between those clocks, 1.84e-8 of the larger within-chart shift. This confirms the near-cancellation through an additional directly evolved clock. It does not bound either individual's count by the precision of their ratio, and the contrast's step sensitivity did not decrease monotonically.

For pressure the corresponding absolute within-chart fractional shifts are +8.17095e-8 and +7.94636e-8. The matched mixed300/G9 contrast changes +2.24584e-9, about 2.75% of the larger within-chart change. Its 4096-step value is −3.98203392e-5 (−0.00398203%). Neither result measures arbitrary changes in the restart, a continuum bound, or a physical density/transfer result.

The 4096-step mixed300 payload retains the original continuity, paired-number and Pauli acceptance and includes the exact checkpoint, source and geometry fingerprints. Replay from the completed payload does not take a new step. The original 1024/2048 joint comparison also reproduces exactly.

```bash
cd kinetics
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
d=2026-11-28; p=results/early-pole-fine-charge2-scalar4-mixed300warm; w=results/early-pole-warm
tmp=$(mktemp); trap 'rm -f "$tmp"' EXIT
gzip -dc "$w-mixed300-x04-4096-checkpoint-$d.json.gz" > "$tmp"
python -W error early_pole_endpoint_warm_x04_mixed300_third_clock.py "$p-x02-state-$d.json.gz" \
 "$p-x04-state-$d.json.gz" "$w-mixed300-actual-prefix-x04-$d.json" \
 "$w-checkpoint-clock-comparison-x04-$d.json" "$w-single-g9-upgrade-$d.json" "$tmp"
python -W error early_pole_endpoint_warm_x04_third_clock_paired_contrast.py \
 "$w-gauss9-x04-third-clock-$d.json" "$w-mixed300-x04-third-clock-$d.json" \
 "$w-checkpoint-clock-comparison-x04-$d.json"
```
