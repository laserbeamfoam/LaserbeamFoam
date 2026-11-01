# twoShots: laserbeamFoam with a long laser pause

A Ti-6Al-4V plate (200 × 100 × 200 µm, 5 µm cells), with gas above, receives two
laser shots of 150 W and 40 µs, at two points 60 µm apart, separated by a pause
of almost 1 ms.

The `longPause` sub-dictionary in `constant/LaserProperties` enables a cheaper
mode for the long laser pauses: once all the metal is solid, only the energy
equation is solved, with the interface frozen and the fluid at rest, with time
steps of up to `maxDeltaT` instead of those of the flow. The solver enters this
thermal mode when all the lasers are off, none is turned on for at least
`minPause`, and all the metal is solid, below the solidus by `solidusMargin`,
and at rest. It goes back to the full equations `leadTime` before the next
laser is turned on. Without `longPause`, laserbeamFoam runs as usual.

In this case:

1. thermal mode until 5 µs before shot 1 (the laser is off for the first
   100 µs);
2. full equations: shot 1 and its solidification;
3. thermal mode: the pause, until 5 µs before shot 2;
4. full equations: shot 2 and its solidification;
5. thermal mode: to the end time (1.4 ms).

Run with `./Allrun` (serial) or `./Allrun parallel` (8 processors,
`system/decomposeParDict`).

## Notes

- In the thermal mode the gas is at rest too, so heat is transferred in the gas
  by conduction only.
- Not available with dynamic mesh refinement or the MHD model. Tested on this
  case only, with one laser.
