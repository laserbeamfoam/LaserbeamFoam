# test_data — End-to-End Post-Processing Test Fixtures

This directory contains the data required by `test_environment.py` to verify
the full post-processing pipeline before running a real calibration.

## Why these files are here

`test_environment.py` does more than check Python imports — it runs
`characterise_meltpool.py` (via `pvpython`) end-to-end on a real simulation
result and compares the extracted meltpool geometry against a known baseline.
This catches broken ParaView installs, wrong `pvpython` versions, and
post-processing regressions that a simple import check would miss.

To do that, it needs actual OpenFOAM field output. The two large files serve
that purpose:

| File | Size | Role |
|---|---|---|
| `0.0007/alpha.metal` | ~210 k lines | Volume-fraction field from a real 316L single-track simulation at t = 0.7 ms; used as input to `pvpython` for VTK extraction |
| `0.0007/meltHistory` | ~210 k lines | Melt-state history field from the same snapshot; also read by `characterise_meltpool.py` |
| `expected_metrics_summary.csv` | 2 lines | Baseline meltpool width/depth/area values that the test compares against |

## Why the files are large

Both fields store one value per mesh cell. The tutorial mesh has ~210 000
cells (300 × 800 × 300 µm domain at ~3 µm cell size), so each field file
has ~210 000 data lines. This cannot be reduced without re-running the
simulation on a coarser mesh and regenerating the baseline CSV, which would
risk making the test less representative of real production runs.

## Regenerating the test data

If the mesh or physics settings change, re-run the tutorial simulation once
and copy the final time-step directory here:

```bash
cp -r runs/<run_id>/<power>W/<timestamp>/  test_data/0.0007/
python characterise_meltpool.py   # run inside test_data/ to get new baseline
cp test_data/metrics_summary.csv test_data/expected_metrics_summary.csv
```

Then commit the updated files.
