# `Plate2D` regression test

This tutorial includes a lightweight regression script:

```bash
./regressionTest.sh
```

The regression test is intended to provide a quick check of the `laserbeamFoam`
thermal and ray-tracing path without running the full tutorial duration.

## What the script does

- copies the tutorial into `regressionTests/main`
- reduces the case `endTime` from `0.01` to `0.001`
- appends a `fieldMinMax` function object for `T`
- runs `./Allrun` in the copied case
- checks a few final values against expected ranges

## Current checks

- final maximum temperature, extracted from
  `postProcessing/TMinMax/0/fieldMinMax.dat`
- final deposited heat reported for `laser0` in `log.laserbeamFoam`
- final deposited heat reported for `laser1` in `log.laserbeamFoam`

These checks are intentionally simple and cheap: they verify that the shortened
case still produces the expected heating response and ray-tracing deposition.

## Notes

- The original tutorial files are not modified during the regression run; all
  changes are made in the copied case under `regressionTests/`.
- You can run the script with `--check-only` to skip rerunning the case and
  only evaluate previously generated regression outputs.
