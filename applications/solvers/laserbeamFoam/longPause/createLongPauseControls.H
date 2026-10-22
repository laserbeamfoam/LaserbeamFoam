// Optional thermal mode for long laser pauses, enabled by a longPause
// sub-dictionary in constant/LaserProperties. During a pause, once all the
// metal is solid and at rest, only the energy equation is solved, with the
// interface frozen and the fluid at rest, until leadTime before the next
// laser is turned on (see checkThermalEntry.H and thermalMode.H)

const bool longPauseActive(laser.found("longPause"));

const dictionary& longPauseDict = laser.subOrEmptyDict("longPause");

// Entry into the thermal mode: all lasers off, no laser turned on for at
// least minPause, no liquid metal, metal below the solidus by solidusMargin
// and metal at rest
const scalar thermalMinPause
(
    longPauseDict.getOrDefault<scalar>("minPause", 20e-6)
);
const scalar thermalMaxLiquidFraction
(
    longPauseDict.getOrDefault<scalar>("maxLiquidFraction", 1e-6)
);
const scalar thermalSolidusMargin
(
    longPauseDict.getOrDefault<scalar>("solidusMargin", 100)
);
const scalar thermalMaxMetalVelocity
(
    longPauseDict.getOrDefault<scalar>("maxMetalVelocity", 1e-3)
);

// Exit from the thermal mode: leadTime before the next laser is turned on
const scalar thermalLeadTime
(
    longPauseDict.getOrDefault<scalar>("leadTime", 5e-6)
);

// Time step in the thermal mode: it starts from initialDeltaT and is limited
// by the largest temperature change per step, by its growth rate and by
// maxDeltaT
const scalar thermalInitialDeltaT
(
    longPauseDict.getOrDefault<scalar>("initialDeltaT", 1e-7)
);
const scalar thermalMaxDeltaT
(
    longPauseDict.getOrDefault<scalar>("maxDeltaT", 1e-5)
);
const scalar thermalMaxTemperatureChange
(
    longPauseDict.getOrDefault<scalar>("maxTemperatureChange", 5)
);
const scalar thermalMaxDeltaTGrowth
(
    longPauseDict.getOrDefault<scalar>("maxDeltaTGrowth", 1.2)
);

if (longPauseActive)
{
    if (thermalMinPause <= thermalLeadTime)
    {
        FatalErrorInFunction
            << "longPause: minPause (" << thermalMinPause
            << ") must be larger than leadTime (" << thermalLeadTime << ")"
            << exit(FatalError);
    }

    if (mthd.valid())
    {
        FatalErrorInFunction
            << "longPause is not available with the electromagnetic (MTHD) "
            << "model" << exit(FatalError);
    }

    if (mesh.dynamic())
    {
        FatalErrorInFunction
            << "longPause is not available with a dynamic mesh"
            << exit(FatalError);
    }

    Info<< "longPause: thermal mode enabled for laser pauses longer than "
        << thermalMinPause << " s" << endl;
}
