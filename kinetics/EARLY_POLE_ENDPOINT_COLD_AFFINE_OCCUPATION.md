# A smooth 0.65% held-occupation ramp reproduces both cold sign reversals

The [constant-occupation audit](EARLY_POLE_ENDPOINT_COLD_NEUTRAL_OCCUPATION_CORRIDORS.md) rules out simultaneous D2 loss and D4 gain on both the full q=[1e−6,1.5e−5] band and its exact upper q=[1e−5,1.5e−5] subband. The [two-zone audit](EARLY_POLE_ENDPOINT_COLD_TWO_ZONE_OCCUPATION.md) shows that a small step can do it. A step is **not** necessary.

Use a common affine profile n(q)=a+b t where t=(q−q_min)/(q_max−q_min) within the retained band. Each native panel has a known integrated positive slope but unknown distribution within the panel. Bound its slope-weighted t by putting that slope at either endpoint. The normalized net in any mask is then a+b E_D[t]−B/D, with both endpoint estimates taken over the exact accepted panels. For D2 loss use the upper mean bound; for D4 gain use the lower mean bound. No within-panel interpolation is needed.

An increasing affine profile with endpoints **2.9046111359e−9 and 2.9234929366e−9** certifies the four signs in both meshes of each chart, even allowing independently unknown positive slope profiles inside every native panel. Their difference is 1.88818007e−11, or **0.6501%** of the lower occupation. The infimum under this support-only affine guarantee is 1.88808007e−11: at the boundary the most restrictive D2 full loss and D4 upper gain constraints meet, while the other two signs have slack. That infimum is about 1.645 times the two-zone infimum. These are different optimization classes, not a contradiction.

The construction is a sensitivity test and uses a small strict sign margin. It neither identifies a live photon distribution nor provides a spatial error bound or physical sign inference.

```bash
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error -m unittest test_early_pole_endpoint_cold_affine_occupation
d=2026-11-28
p=results/early-pole-cold
python -W error early_pole_endpoint_cold_affine_occupation.py results \
  "$p-neutral-occupation-corridors-$d.json" "$p-two-zone-occupation-$d.json" |
  cmp - "$p-affine-occupation-$d.json"
```
