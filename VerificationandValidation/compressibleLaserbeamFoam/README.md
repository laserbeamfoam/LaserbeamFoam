# Verification suite: set-up fix (10 Nov)

## What went wrong in 11–16

The 3D cases keep their base-mesh fields in `0.orig`. Each case's `Allrun` copies it to `0/`, lets the solver refine for two steps, then rewrites the fields on the refined mesh.

In your copy `0.orig` was missing, so the copy failed and everything after it failed too. The generators do write it, so it went missing on the way to the machine that ran the suite. Two common causes:

- **A `*.orig` ignore rule.** Many global gitignores have one, for merge backups, and so do some rsync or backup exclude lists. It drops exactly `0.orig` and leaves the 1D cases' `0/` alone, which matches what you saw.
- **The LaserbeamFoam repo's `.gitignore`.** It ignores `0*` … `9*`. That matches `0.orig` and also every case folder (`01_…` to `19_…`). If the suite lives in that checkout, commit only the generators and `tools/`, and regenerate the cases.

## What changed

| file | change |
|---|---|
| `gen_cases3d.py` | New `--orig CASE` mode rewrites `0.orig` from the case's spec (byte-identical to what the generator writes). The case `Allrun` calls it when `0.orig` is missing. The `Allrun` also stops at the first failed step and names its log (refining steps, initialisation, decomposePar), and exits with the solver's status. |
| `gen_cases.py` | The 1D `Allrun` exits with the solver's status, so a crash counts as a failure. |
| `Allrun` (top level) | A failed case prints the last 15 lines of its newest log. |

The case physics, meshes, fields and numerics are unchanged.

## To use

```
cp gen_cases.py gen_cases3d.py Allrun  <your clbf_verification>/
cd <your clbf_verification>
python3 gen_cases.py cases; python3 gen_cases3d.py cases   # the case Allruns come from the generators
./Allrun 11 12 13 14 15 16 19
```

Regenerating also writes `0.orig` again, so 11–16 run straight away.

## Checked here (v1912)

- With `0.orig` deleted, 13 and 14 rebuilt it and ran to the end through the top-level `Allrun`.
- With `0/T` deleted, case 05 was reported as failed, followed by its FatalError.
