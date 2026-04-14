//! put the usual header here:
/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of solids4foam.

    solids4foam is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    solids4foam is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with solids4foam.  If not, see <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "beamShaper.H"
#include "volFields.H"
#include "dictionary.H"
#include "interpolationTable.H"
#include "OFstream.H"
#include "IFstream.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// defineTypeNameAndDebug(BeamShape, 0);

// * * * * * * * * * * * * Private Member Functions  * * * ** * * * * * * * * * //

BeamShape::BeamShape
(
    const fileName& shapeFile,
    const fileName& scheduleFile,
    const fvMesh& mesh,
    const dictionary& laserProperties // This line accepts the LaserProperties dictionary
)
: 
    shapeFile_(shapeFile), // Initialize shapeFile_
    scheduleFile_(scheduleFile),
    mesh_(mesh),           // Initialize mesh_ if it is needed later
    laserProperties_(laserProperties) // Properly initialize laserProperties_ here
{
    // read the laserProperties file for the debug flag
    debug = laserProperties.lookupOrDefault<bool>("debug", false);
    
    // read the LaserProperties dictionary here for the pixel size
    
    pixelSize_ = laserProperties.lookupOrDefault<scalar>("pixelSize", 2.5814851E-6); // default for 1.5 m free space propagation
    pixelArea_ = pow(pixelSize_, 2);
    // set the imageLength_
    imageLength_ = dim_ * pixelSize_;
    
    addressShape_ = List<scalar>(pow(dim_,2), 0.0); // Assign storage for the current address shape
    beamShape_ = List<scalar>(pow(dim_,2), 0.0);   // Use Foam::zero to initialize the field
    beamBox_   = List<scalar>(pow(dim_exp,2), 0.0);  // Assign storage for the expanded field
    GaussianField_ = List<scalar>(pow(dim_,2), 0.0); // Assign storage for the Gaussian field
    // beamShapeAddresses_ = List<Pair<scalar>>;
    beamShapeGenerator();

    IOdictionary beamShapeDict
    (
        IOobject
        (
            scheduleFile_,
            mesh.time().constant(),  // Assuming "constant" time directory
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    // Step 1: Read the list of shape names from the namesDict file
    beamShapesList_ = readShapeNames(beamShapeDict);
    dummyAddress_ = Pair<label>(1936588847, 1936588847);
    if (debug)
    {
        Info << "Number of pair lists to read: " << beamShapesList_.size()
            << endl;
    }

    // Step 2: Loop over the list of names and read each pair list from addressListNamesDict
    forAll(beamShapesList_, i)
    {
        const word& listName = beamShapesList_[i];
        if (beamShapeDict.found(listName))
        {
            List<Pair<label>> addressList = readAddressesByName(beamShapeDict, listName);

            if (debug)
            {
                Info << "Reading pair list: " << listName << " with "
                    << addressList.size() << " entries" << endl;
            }
            beamShapeLibrary_.append(addressList);
        }
        else
        {
            FatalErrorInFunction << "Pair list " << listName << " not found in dictionary" << exit(FatalError);
        }
    }

    beamShapeScheduler(beamShapeDict);

}

BeamShape::~BeamShape()
{
    // Destructor implementation (if needed)
}

List<word> BeamShape::readShapeNames
(
    const dictionary& dict
)
{
    return dict.lookup("shapeNames");
}

// Function to read a list of Pair<label> from a dictionary by name
List<Pair<label>> BeamShape::readAddressesByName
(
    const dictionary& dict, 
    const word& listName
)
{
    return dict.lookup(listName);
}

void BeamShape::printBeamShape
(
    const fileName& filename, 
    List<scalar>& arrayData, 
    label dim
)
{
    // Open the file
    OFstream file(filename);

    // Loop over rows and columns
    for (label row = 0; row < dim; ++row)
    {
        for (label col = 0; col < dim; ++col)
        {
            const label i = row * dim + col;
            scalar value = getValue(arrayData, i);

            file << value;

            if (col < dim - 1)
            {
                // Not the last column, add a comma
                file << ",";
            }
        }
        // End of row, add a newline
        file << "\n";
    }

    // Close the file (optional, destructor will handle it)
    // file.close();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void BeamShape::beamShapeGenerator()
{
    
    if (debug)
    {
        Info << "Reading primitive shape from: " << shapeFile_ << endl;
    }

    // Open the CSV file
    IFstream primitiveFile(shapeFile_);

    if (!primitiveFile.good())
    {
        FatalErrorInFunction << "Cannot open primitive shape file: " << shapeFile_ << exit(FatalError);
    }

    // Read data into primitiveShape_
    DynamicList<scalar> dataList;

    while (!primitiveFile.eof())
    {   
        scalar value;
        primitiveFile >> value;
        dataList.append(value);
        // std::cin.get();
    }
    addressShape_.transfer(dataList);

    // forAll(addressShape_, i)
    // {
    //     if (addressShape_[i] < 0.5) addressShape_[i] = 0.0;
    // }

    if (debug)
    {
        printBeamShape("addressShape.csv", addressShape_, dim_);
    }

    normalizeField(addressShape_);

    generateGaussianFilter(sigma);

    applyGaussianFilter(addressShape_, 1);

    generateExpandedBeamShape(7);

    //     List<scalarField> beamData_collect(Pstream::nProcs());
    //     beamData_collect[Pstream::myProcNo()] = beamData_;
    //     Pstream::gatherList(beamData_collect); // Gather the data onto the master processor
    //     Pstream::scatterList(beamData_collect);  // Scatter the data back to all processors
    if (debug)
    {
        Info << "Finished beamShapeGenerator." << endl;
    }

}

void BeamShape::generateExpandedBeamShape
(
    const label num_copies_x
)
{
    if (debug)
    {
        Info << "Generating expanded beam shape." << endl;
    }
    List<scalar> templateBeamShape = List<scalar>(41209, 0.0);
    // getting the central mappable image from the beam shape
    
    forAll (templateBeamShape, i)
    {
        // mapping to the central 203 x 203 pixels of the beamShape_
        label x_corner = 175;
        label z_corner = 175;
        Pair<label> coords_i = coords(i, dim_inp);
        label x_offset = coords_i.first()+x_corner;
        label z_offset = coords_i.second()+z_corner;
        // scalar value = ;
        setValue(templateBeamShape, i, getValue(addressShape_, index(Pair<label>(x_offset, z_offset), dim_)));
    }

    if (debug)
    {
        printBeamShape("templateBeamShape.csv", templateBeamShape, dim_inp);
    }

    // Close the file (optional, destructor will handle it)

    // take the central 203 * 203 pixels from the beam shape
    // and repeat it 5 times in the x and y directions
    // with a 1 pixel gap between each copy and starting 1 pixel in from the edge
    // to avoid boundary effects

    // TODO: convert all references to x and y in the beamShape code refer to x and z ()

    for (label row = 0; row < num_copies_x; ++row)
    {
        for (label col = 0; col < num_copies_x; ++col)
        {
            forAll(templateBeamShape, i)
            {
                Pair<label> coords_i = coords(i, dim_inp);
                label x_offset = 1 + col * (dim_inp+1);
                label z_offset = 1 + row * (dim_inp+1);
                // set the value
                setValue(beamBox_, index(
                    Pair<label>(coords_i.first() + x_offset, coords_i.second() + z_offset), dim_exp), 
                    getValue(templateBeamShape, i));
            }
        }
    }
    forAll(beamBox_, i)
    {
        if (beamBox_[i] < SMALL) beamBox_[i] = SMALL;
    }

    if (debug)
    {
        printBeamShape("expandedBeamShape.csv", beamBox_, dim_exp);
    }
    
}

void BeamShape::beamShapeScheduler
(
    const dictionary& beamShapeDict
)
{
    shapeFrequency_ = List<scalar>(beamShapeDict.lookup("shapeFrequency"));
    shapeDuration_ = List<scalar>(beamShapeDict.lookup("shapeDuration"));
    cycleDuration_ = sum(shapeDuration_);
    if (debug)
    {
        Info << "Shape duration: " << shapeDuration_ << endl;
        Info << "Cycle duration: " << cycleDuration_ << endl;
    }

    label addressCount = 0;
    // Clear previous scheduling
    addressDuration_.clear();
    beamShapeIndex_.clear();
    beamShapeTiming_.clear();
    beamShapeAddresses_.clear();

    forAll(beamShapeLibrary_, i)
    {
        const List<Pair<label>>& addressList = beamShapeLibrary_[i];
        label n_addresses = addressList.size();
        scalar freq = shapeFrequency_[i];
        scalar shapeDur = shapeDuration_[i];

        if (n_addresses == 1)
        {
            // Only one address in this shape
            scalar addrDur = shapeDur;
            addressDuration_.append(addrDur);
            addressCount++;
            beamShapeIndex_.append(addressCount);
            beamShapeTiming_.append(addrDur);
            beamShapeAddresses_.append(addressList[0]);
            // No dummy address or interchange
        }
        else
        {
            // Multiple addresses: calculate address duration and repeat
            scalar addrDur = 1.0 / (n_addresses * freq) - addressInterchange_;
            addressDuration_.append(addrDur);
            // How many full cycles of all addresses fit in the shape duration?
            scalar cycleTime = n_addresses * (addrDur + addressInterchange_);
            label n_cycles = std::floor(shapeDur / cycleTime + SMALL); // add SMALL for floating point safety
            scalar usedTime = n_cycles * cycleTime;
            scalar remainingTime = shapeDur - usedTime;
            // Schedule full cycles
            for (label c = 0; c < n_cycles; ++c)
            {
                forAll(addressList, j)
                {
                    addressCount++;
                    beamShapeIndex_.append(addressCount);
                    beamShapeTiming_.append(addrDur);
                    beamShapeAddresses_.append(addressList[j]);
                    beamShapeTiming_.append(addressInterchange_);
                    beamShapeAddresses_.append(dummyAddress_);
                }
            }
            // Handle any remaining time (partial cycle)
            label j = 0;
            while (remainingTime > 0.0 && j < n_addresses)
            {
                scalar t_this = std::min(addrDur, remainingTime);
                addressCount++;
                beamShapeIndex_.append(addressCount);
                beamShapeTiming_.append(t_this);
                beamShapeAddresses_.append(addressList[j]);
                remainingTime -= t_this;
                if (remainingTime > 0)
                {
                    scalar t_inter = std::min(addressInterchange_, remainingTime);
                    beamShapeTiming_.append(t_inter);
                    beamShapeAddresses_.append(dummyAddress_);
                    remainingTime -= t_inter;
                }
                ++j;
            }
        }
    }

    // pick up numerical errors in the cycle accumulation
    scalar timeSum = 0.0;
    forAll(beamShapeTiming_, i)
    {
        timeSum += beamShapeTiming_[i];
        if (debug)
        {
            Info << "Beam shape timing at " << i << ": Address: "
                << beamShapeAddresses_[i] << " timing: "
                << beamShapeTiming_[i] << " sum: " << timeSum << endl;
        }
    }
    if (mag(timeSum - cycleDuration_) > SMALL) FatalErrorInFunction << "Beam shape timing =/= cycleDuration_." << exit(FatalError);

    if (debug)
    {
        Info << "Finished scheduling beam shapes." << endl;
    }
}

void BeamShape::getBeamShape
(
    scalar currentTime,
    scalar timeStep
)
{
    /***********************************\
     * this function will be called by *
     * laserHeatSource::updateDeposition* 
    \***********************************/

    shapeOverInterval(currentTime, currentTime + timeStep);

}

label BeamShape::currentAddressIndex
(
    scalar currentInCycleTime
)
{
    /** this function takes a time, and returns the 
     * index of the current beam shape address in 
     * the beamShapeAddresses_ list based on the 
     * beamShapeTiming_ list
     **/

    // Info << "Getting current address index." << endl;
    // scalar inCycleTime = fmod(currentTime, cycleDuration_);
    // Info << "Current time: " << currentTime << " inCycleTime: " << inCycleTime << endl;
    // find the index of the current time in the beamShapeTiming_ list
    label index = 0;
    scalar timeSum = 0.0;
    
    forAll(beamShapeTiming_, i)
    {
        timeSum += beamShapeTiming_[i];
        // Info << "Beam shape timing at " << i << ": " << beamShapeTiming_[i] << " sum: "<< timeSum << endl;
        if (timeSum > currentInCycleTime)
        {
            // Info << "Found index: " << i << endl;
            index = i;
            break;
        }
    }

    // Info << "Current time: " << currentTime << " inCycleTime: " << inCycleTime << endl;

    return index;
}

Pair<label> BeamShape::currentBeamShapeAddress
(
    scalar currentTime
)
{
    return beamShapeAddresses_[currentAddressIndex(currentTime)];
}

/***********************************************\
|* BeamShape computed over a time interval     *|
\***********************************************/
void BeamShape::shapeOverInterval(scalar currentTime, scalar nextTime)
{
    const scalar tolerance = max(SMALL, scalar(1e-15));
    const scalar timeStep = nextTime - currentTime;

    List<label> addressIndices;
    List<scalar> addressWeights;

    if (timeStep <= tolerance)
    {
        scalar inCycleTime = fmod(currentTime, cycleDuration_);
        if (inCycleTime < 0) inCycleTime += cycleDuration_;
        if (cycleDuration_ - inCycleTime < tolerance)
        {
            inCycleTime = 0.0;
        }

        label currentIndex = currentAddressIndex(inCycleTime);
        addressIndices.append(currentIndex);
        addressWeights.append(1.0);
        generateMeanImage(addressIndices, addressWeights);
        return;
    }

    List<scalar> cumulativeTimes(beamShapeTiming_.size() + 1, scalar(0));
    forAll(beamShapeTiming_, i)
    {
        cumulativeTimes[i + 1] = cumulativeTimes[i] + beamShapeTiming_[i];
    }

    scalar remaining = timeStep;
    scalar absoluteTime = currentTime;

    while (remaining > tolerance)
    {
        scalar inCycleTime = fmod(absoluteTime, cycleDuration_);
        if (inCycleTime < 0) inCycleTime += cycleDuration_;
        if (cycleDuration_ - inCycleTime < tolerance)
        {
            inCycleTime = 0.0;
        }

        label currentIndex = currentAddressIndex(inCycleTime);
        scalar addressEnd   = cumulativeTimes[currentIndex + 1];
        scalar addressStart = addressEnd - beamShapeTiming_[currentIndex];
        scalar timeLeftInAddress = addressEnd - inCycleTime;
        scalar safeTimeLeft = max(timeLeftInAddress, tolerance);
        scalar contribution = min(safeTimeLeft, remaining);

        if (contribution <= tolerance)
        {
            absoluteTime += safeTimeLeft;
            remaining = max(remaining - safeTimeLeft, scalar(0));
            continue;
        }

        const Pair<label>& currentAddress = beamShapeAddresses_[currentIndex];

        if (currentAddress == dummyAddress_)
        {
            // Redistribute dummy contributions to neighbouring real addresses
            const label nAddresses = beamShapeAddresses_.size();
            label prevIndex = currentIndex;
            label nextIndex = currentIndex;
            bool foundPrev = false;
            bool foundNext = false;

            for (label step = 0; step < nAddresses; ++step)
            {
                prevIndex = (prevIndex == 0 ? nAddresses - 1 : prevIndex - 1);
                if (beamShapeAddresses_[prevIndex] != dummyAddress_)
                {
                    foundPrev = true;
                    break;
                }
            }

            for (label step = 0; step < nAddresses; ++step)
            {
                nextIndex = (nextIndex + 1) % nAddresses;
                if (beamShapeAddresses_[nextIndex] != dummyAddress_)
                {
                    foundNext = true;
                    break;
                }
            }

            if (!foundPrev || !foundNext)
            {
                FatalErrorInFunction
                    << "Unable to locate neighbouring beam-shape addresses"
                    << " for dummy slot at index " << currentIndex << exit(FatalError);
            }

            scalar localStart = inCycleTime;
            scalar localEnd = inCycleTime + contribution;
            const scalar intervalLength = max(addressEnd - addressStart, tolerance);
            scalar averageFraction = ((localStart + localEnd)*0.5 - addressStart)/intervalLength;
            averageFraction = min(max(averageFraction, scalar(0)), scalar(1));

            scalar weightPrev = (1.0 - averageFraction)*contribution;
            scalar weightNext = averageFraction*contribution;

            if (weightPrev > tolerance)
            {
                addressIndices.append(prevIndex);
                addressWeights.append(weightPrev);
            }

            if (weightNext > tolerance)
            {
                addressIndices.append(nextIndex);
                addressWeights.append(weightNext);
            }
        }
        else
        {
            addressIndices.append(currentIndex);
            addressWeights.append(contribution);
        }

        remaining    -= contribution;
        absoluteTime += contribution;
    }

    List<scalar> aggregatedWeights(beamShapeTiming_.size(), scalar(0));
    forAll(addressIndices, i)
    {
        aggregatedWeights[addressIndices[i]] += addressWeights[i];
    }

    List<label> finalIndices;
    List<scalar> finalWeights;
    forAll(aggregatedWeights, idx)
    {
        if (aggregatedWeights[idx] >= tolerance)
        {
            finalIndices.append(idx);
            finalWeights.append(aggregatedWeights[idx]);
        }
    }

    if (finalIndices.empty())
    {
        FatalErrorInFunction
            << "Failed to accumulate beam shape contributions for interval ["
            << currentTime << ", " << nextTime << "] with time step "
            << timeStep << "." << exit(FatalError);
    }

    if (debug) Info << "Address weights: " << finalWeights << endl;

    generateMeanImage(finalIndices, finalWeights);

    if (debug)
    {
        fileName filename = "beamShape_" + std::to_string(currentTime) + "_" + std::to_string(nextTime) + ".csv";
        printBeamShape(filename, beamShape_, dim_);
    }
}

void BeamShape::generateMeanImage(
    List<label> addressIndices,
    List<scalar> addressWeights
)
{
    // Set all the values of the beamShape_ to zero
    forAll(beamShape_, i)
    {
        setValue(beamShape_, i, 0.0);
    }

    // 33x33 2D array for weights, mapped from address pairs (-16..16)
    const label gridMin = -16;
    const label gridMax = 16;
    const label gridSize = gridMax - gridMin + 1; // 33
    List<List<scalar>> addressWeight2D(gridSize);
    for (label i = 0; i < gridSize; ++i)
        addressWeight2D[i] = List<scalar>(gridSize, 0.0);

    // Compile weights into 2D array
    forAll(addressIndices, idx)
    {
        Pair<label> address = beamShapeAddresses_[addressIndices[idx]];
        label x = address.first();
        label y = address.second();
        // Map (-16..16) to (0..32)
        label xi = x - gridMin;
        label yi = y - gridMin;
        if (xi >= 0 && xi < gridSize && yi >= 0 && yi < gridSize)
            addressWeight2D[xi][yi] += addressWeights[idx];
    }

    // For each address with nonzero weight, call beamShapeFromAddress once
    scalar effectiveWeight = 0.0;
    for (label xi = 0; xi < gridSize; ++xi)
    {
        for (label yi = 0; yi < gridSize; ++yi)
        {
            scalar weight = addressWeight2D[xi][yi];
            if (weight >= SMALL)
            {
                Pair<label> address(xi + gridMin, yi + gridMin);
                beamShapeFromAddress(address);
                effectiveWeight += weight;
                forAll(beamShape_, j)
                {
                    setValue(beamShape_, j, getValue(beamShape_, j) + weight * getValue(addressShape_, j));
                }
            }
        }
    }

    // Normalize by total weight
    if (effectiveWeight > VSMALL)
    {
        forAll(beamShape_, i)
        {
            setValue(beamShape_, i, getValue(beamShape_, i) / effectiveWeight);
        }

    normalizeField(beamShape_); // make integral unity
    }
    else
    {
        FatalErrorInFunction
            << "Total weight for mean image is zero or near zero." << nl
            << "Indices: " << addressIndices << nl
            << "Weights: " << addressWeights << exit(FatalError);
    }
}

void BeamShape::beamShapeFromAddress
(
    const Pair<label>& beamShapeAddress
)
{
    /*****************************************************\
    This function take a beamShapeAddress (pair of labels)
    and produces a 550x550 beam shape from the expanded, 
    setting the current beam shape to the result.
    Adapt from this python code, using the c++ structures:
    corner_x = (self.centre_x_p - 275 + 103) - int(106*2*(peak_address[0])/self.grid_size)
    corner_y = (self.centre_y_p - 275 + 103) - int(106*2*(peak_address[1])/self.grid_size)
    # Ensure the indices are within bounds
    if (corner_y + 550 <= self.expanded.shape[0]) and (corner_x + 550 <= self.expanded.shape[1]):
        temporary = self.expanded[corner_y:corner_y+550, corner_x:corner_x+550]
        # print('Temporary shape:', temporary.shape)
        # temporary = self.apply_gaussian(temporary, self.sigma)  # Ensure this is 2D
    else:
        print('Error: Indices out of bounds')
        return np.zeros((550, 550))
    \*****************************************************/
    // set to zero 
    forAll(addressShape_, i)
    {
        setValue(addressShape_, i, 0.0);
    }
    // check is the address is the dummy address
    if (beamShapeAddress != dummyAddress_)
    {
        
        // Info << "Generating beam shape from address: " << beamShapeAddress << endl;
        // Info << "Centre of the beam shape: " << beamShapeAddress.first() << ", " << beamShapeAddress.second() << endl;
                       // image centroid - beamShape/2     
        label corner_row =
                ( dim_exp / 2 ) + ( beamShapeAddress.second() * dim_inp / dim_grd )
            - ( dim_ / 2 );
        label corner_col =
                ( dim_exp / 2 ) - ( beamShapeAddress.first() * dim_inp / dim_grd )
            - ( dim_ / 2 );

        if (debug) 
        {
            Info << "Address: " << beamShapeAddress << endl;
            Info << "corner_row: " << corner_row
                 << ", corner_col: " << corner_col << endl;
        }
        forAll(addressShape_, i)
        {
            Pair<label> coords_i = coords(i, dim_);
            label beamBox_row = corner_row + coords_i.first();
            label beamBox_col = corner_col + coords_i.second();
            setValue
            (
                addressShape_,
                i,
                getValue
                (
                    beamBox_,
                    index(Pair<label>(beamBox_row, beamBox_col), dim_exp)
                )
            );
        }

        // Normalize the beam shape
        normalizeField(addressShape_);
        applyGaussianFilter(addressShape_, false);
        normalizeField(addressShape_);

        if (debug)
        {
            fileName filename =
                "beamShape_" + std::to_string(beamShapeAddress.first())
              + "_" + std::to_string(beamShapeAddress.second()) + ".csv";
            printBeamShape(filename, addressShape_, dim_);
        }
    }
    else
    {
        FatalErrorInFunction
            << "Encountered dummy beam-shape address during generation"
            << exit(FatalError);
    }
    
}

void BeamShape::generateGaussianFilter
(
    const scalar sigma
)
{
    if (debug)
    {
        Info << "Computing the Gaussian Field." << endl;
    }
    // Generate a Gaussian filter to apply to the beam shape
    forAll(GaussianField_, i)
    {
        // Info << "Computing Gaussian Field at: " << i << endl;
        Pair<label> coords_i = coords(i, 550);
        scalar x_i = coords_i.first()-275;
        scalar z_i = coords_i.second()-275;
        scalar g_i = Foam::exp(-(x_i*x_i + z_i*z_i)/(2*sigma*sigma))
                    / (2*Foam::constant::mathematical::pi*sigma*sigma);
        setValue(GaussianField_, i, g_i);
    }
    // Info << "Done." << endl;
    // Normalize the Gaussian (unit volume)
    if (debug)
    {
        Info << "Normalizing the Gaussian." << endl;
    }
    normalizeField(GaussianField_);
    if (debug)
    {
        printBeamShape("GaussianField.csv", GaussianField_, dim_);
    }
    // Info << "Done." << endl;
}

void BeamShape::applyGaussianFilter
(
    List<scalar>& arrayData,
    bool inverse
) 
{   
    // Info << "Applying Gaussian filter. Inverse = " << inverse << endl;
    // Info << "Applying field." << endl;
    forAll(arrayData, i)
    {
        scalar value = getValue(arrayData, i);
        scalar g = getValue(GaussianField_, i);
        // Info << "At i = " << i << ", Value: " << value << ", Gaussian: " << g << endl;
        if (inverse) 
        {
            scalar temp = value/g;
            setValue(arrayData, i, temp);
        }
        else 
        {
            scalar temp = value*g;            
            setValue(arrayData, i, temp);
        }
    }
    // Info << "Done." << endl;
}

void BeamShape::normalizeField(List<scalar>& arrayData) {
    // Calculate the sum of all values in arrayData
    scalar sum = 0.0;
    forAll(arrayData, i) {
        sum += getValue(arrayData, i);
    }
    // Info << "Sum of elements: " << sum << endl;

    // Check if sum is close to zero to avoid division errors
    if (mag(sum) < SMALL) {
        FatalErrorInFunction << "Normalization failed: sum of elements is zero or near zero." << exit(FatalError);
    }

    // Normalize each element to make the integral (sum) equal to 1
    forAll(arrayData, i) {
        setValue(arrayData, i, getValue(arrayData, i) / sum);
    }

}
inline label BeamShape::index
(
    const Pair<label>& coords,
    const label dim
)
{
    return coords.first()*dim + coords.second();
}

Pair<label> BeamShape::coords
(
    const label index,
    const label dim
)
{
    return Pair<label>(index/dim, index%dim);
}

void BeamShape::setValue(
    List<scalar>& arrayData,
    label dim,
    const Pair<label>& coords,
    const scalar value
)
{
    arrayData[index(coords, dim)] = value;
}

void BeamShape::setValue(
    List<scalar>& arrayData,
    const label index,
    scalar value
)
{
    arrayData[index] = value;
}

scalar BeamShape::getValue(
    List<scalar>& arrayData,
    const Pair<label>& coords, // Make this parameter a const reference
    label dim
) 
{
    return arrayData[index(coords, dim)];
}

scalar BeamShape::getValue
(
    List<scalar>& arrayData,
    label index
) 
{
    return arrayData[index];
}

scalar BeamShape::interpolateValue
(
    scalar x_target,
    scalar z_target,
    scalar localCellSize,
    label dim
) 
{
    // for our purpose, multilinear interpolation is sufficient
    // we will use the four nearest neighbours to the target point
    // and interpolate the value based on the distance to each point
    
    // the scale of the beam shape is 550*2.5814851E-6
    // the target point is in the range of -275*2.5814851E-6 to 275*2.5814851E-6

    if (localCellSize < SMALL)
    {
        if (debug)
        {
            Info << "Warning: Very small localCellSize: " << localCellSize 
                 << " at point (" << x_target << ", " << z_target << ")" << endl;
        }
        // Use a minimum safe cell size to prevent numerical issues
        localCellSize = SMALL;
    }

    // Calculate the areal scaling factor based on the local cell size
    // Define the bounds for the beam shape (original pixel size range)
    scalar halfBeamRange = 275 * pixelSize_.value();

    // Check if the target point is within the bounds of the beam shape
    if (fabs(x_target) > halfBeamRange || fabs(z_target) > halfBeamRange) 
    {
        return 0.0;
    }


    // Convert the target point into beam shape units (based on pixel size)
        scalar row = 275 + (-z_target / pixelSize_.value());
        scalar col = 275 + (x_target / pixelSize_.value());

    // Determine the four nearest neighbors
        label row1 = std::floor(row);
        label row2 = std::ceil(row);
        label col1 = std::floor(col);
        label col2 = std::ceil(col);

    // Ensure indices are within bounds
        row1 = std::max(0, std::min(row1, dim -1));
        row2 = std::max(0, std::min(row2, dim -1));
        col1 = std::max(0, std::min(col1, dim -1));
        col2 = std::max(0, std::min(col2, dim -1));

    // Debugging output 
    if (debug)
    {
        Info << "Target point: (" << x_target << ", " << z_target << ")" << endl;
           Info << "Beam shape point: (row=" << row << ", col=" << col
               << ")" << endl;
           Info << "Nearest neighbors: (" << row1 << ", " << col1 << "), ("
               << row1 << ", " << col2 << ")" << endl;
           Info << "                   (" << row2 << ", " << col1 << "), ("
               << row2 << ", " << col2 << ")" << endl;
    }   

    // Create pairs for the coordinates of the four nearest neighbors
        Pair<label> coords_v11(row1, col1);
        Pair<label> coords_v12(row1, col2);
        Pair<label> coords_v21(row2, col1);
        Pair<label> coords_v22(row2, col2);

    // Get the values at the four nearest neighbors
    scalar v11 = getValue(beamShape_, coords_v11, dim);
    scalar v12 = getValue(beamShape_, coords_v12, dim);
    scalar v21 = getValue(beamShape_, coords_v21, dim);
    scalar v22 = getValue(beamShape_, coords_v22, dim);

    // Perform bilinear interpolation
    scalar v = (v11 * (row2 - row) * (col2 - col) +
                v21 * (row - row1) * (col2 - col) +
                v12 * (row2 - row) * (col - col1) +
                v22 * (row - row1) * (col - col1));

    //! Catch the (literally) edge case where the interpolated value is negative 
    // This typically occurs where the mesh is very coarse around the edge of the beam shape
    // Scale the interpolated value based on the areal scaling factor
    // return (v * arealScalingFactor / pixelArea_.value()); 
    // (shouldn't happen often, since we refine the mesh around the beam shape)
    // else, return zero
    if (v < 0.0)
    {
        if (fabs(v) > 1e-8)
        {
            Info << "Interpolated value is negative: " << v << endl;
        }
        return 0.0;
    }
    scalar result = v / pixelArea_.value();
    return result;
}

} // End namespace Foam
