Info << "Setting nucleation sites" << endl;

Info<< "Reading Nucleation Parameters from nucleationDict\n" << endl;
IOdictionary nucleationDict
(
    IOobject
    (
        "nucleationDict",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

// Bulk Nucleation Parameters
dimensionedScalar Tu_mean("Tu_mean", nucleationDict);



Info << "Mean undercooling: " << Tu_mean.value() << " K" << endl;  

dimensionedScalar Tu_stdev("Tu_stdev", nucleationDict);




Info << "Standard deviation undercooling: " << Tu_stdev.value() << " K" << endl;

dimensionedScalar n_max("n_max", nucleationDict);




Info << "n max: " << n_max.value() << " sites/m3" << endl;

// Check if mesh is 2D, if so scale volumetric n_max
// 0.7236012546 = sqrt(pi/6)
if(mesh.nGeometricD() == 2){
    Info << "Mesh is 2D. Converting 3D n_max = " << n_max.value() << " sites/m3 to 2D n_max = ";
    n_max.value() = Foam::pow( n_max.value() / 0.7236012546 , (2.0/3.0) ) ;
    Info << n_max.value() << " sites/m2" << endl;
}

scalar nucDistFactor // how many interface thicknesses each nuc site should be from each other
(
    readScalar(nucleationDict.lookup("nucDistFactor"))
);

Info << "Setting nucleation sites with activation undercoolings based on Gaussian distribution" << endl;

// BULK NUCLEATION
// n_max is the maximum density of nucleation sites possible.
scalar meshVolume = gSum(mesh.V());

Pout << "Proc " << Pstream::myProcNo() << ":\tnCells = " << mesh.nCells() << endl;
label totalCellsGlobal(mesh.nCells());
// Need to sum up over all partitions
reduce(totalCellsGlobal,sumOp<label>());

Info << "\tSUM OP nCells: " << totalCellsGlobal << endl;

// Total possible number of nuc sites is meshVolume * n_max
scalar N_maxScalar ( n_max.value() * meshVolume);

Info << "Before inting N_max: " << N_maxScalar 
<< "\tINT_MAX = " << INT_MAX 
<< "\t(" << 100.0*(N_maxScalar/INT_MAX) << " %)" 
<< endl;

if (N_maxScalar > INT_MAX)
{
    FatalError << "Total number of potential nucleation sites ("<<N_maxScalar <<")"
    << " exceeds maximum integer possible (" << INT_MAX 
    << ") [" << 100.0*(N_maxScalar/INT_MAX) << " %]" 
    << abort(FatalError);
}

if (N_maxScalar > totalCellsGlobal)
{
    FatalError << "Total number of potential nucleation sites ("<<N_maxScalar <<")"
    << " exceeds number of cells in the mesh (" << totalCellsGlobal 
    << ") [" << 100.0*(N_maxScalar/totalCellsGlobal) << " %]" 
    << abort(FatalError);
}

int N_max (N_maxScalar);


Info << "After inting N_max: " << N_max << endl;

// Proportion of sites that can nucleate
scalar pSite = scalar(N_max) / totalCellsGlobal;

Info << "Total mesh volume: " << meshVolume 
<< "\tn_max: " << n_max << " sites/m3" 
<< "\tN_max: " << N_max << " sites (" << pSite * 100.0 << "%% of cells)" 
<< endl;


label ni_fieldsNeeded(N_max);

Info<< "Creating nucleation undercooling field (nucTu)" << endl;
volScalarField nucTu
(
    IOobject
    (
        "nucTu",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ, // Could set to READ_IF_PRESENT if desired
        IOobject::AUTO_WRITE // Could be manually written instead
    ),
    mesh,
    dimensionedScalar("nucTu",dimensionSet(0, 0, 0, 1, 0),-100.0) // -100 is a dummy value
);


Info << "Marking nucleation sites (nucSite):\t"
    << "-1 = Not a site\t"
    << "0 = Site that is NOT active\t"
    << "1 = Sit that IS active."
    << endl;

volScalarField nucSite
(
    IOobject
    (
        "nucSite",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE // Could be manually written instead
    ),
    mesh,
    dimensionedScalar("nucSite",dimensionSet(0, 0, 0, 1, 0),-1.0)
);




// Possible solution:
// Generate a list of random locations within the GLOBAL mesh dimensions
//      - Use reduce max and min across processors for X,Y,Z
//      - Generate random location: p = [rand(Xmin,Xmax), rand(Ymin,Ymax)...
//		- Check if this is too close to any existing sites
// 			- if invalid, generate a new one, repeat until valid
//      - Sync list across processors, so each processor is using the same points
//		- Distribute these sites across processors using findCell


// Get global mesh size
boundBox meshBox(mesh.points());

// Get local mesh sizes
scalar Xmin = meshBox.min().component(0);
scalar Ymin = meshBox.min().component(1);
scalar Zmin = meshBox.min().component(2);

scalar Xmax = meshBox.max().component(0);
scalar Ymax = meshBox.max().component(1);
scalar Zmax = meshBox.max().component(2);

// NOT NEEDED, boundBox gets the global mesh size, but it can't hurt
// // use MPI reduce, to find the global mesh sizes
scalar XminG = Xmin; reduce(XminG,minOp<scalar>());
scalar YminG = Ymin; reduce(YminG,minOp<scalar>());
scalar ZminG = Zmin; reduce(ZminG,minOp<scalar>());
scalar XmaxG = Xmax; reduce(XmaxG,maxOp<scalar>());
scalar YmaxG = Ymax; reduce(YmaxG,maxOp<scalar>());
scalar ZmaxG = Zmax; reduce(ZmaxG,maxOp<scalar>());

// 
Pout << "Processor " << Pstream::myProcNo() 
<< endl
<< "\tx:(" << Xmin << "," << Xmax << ")"
<< "\tX:(" << XminG << "," << XmaxG << ")"
<< endl
<< "\ty:(" << Ymin << "," << Ymax << ")"
<< "\tY:(" << YminG << "," << YmaxG << ")"
<< endl
<< "\tz:(" << Zmin << "," << Zmax << ")"
<< "\tZ:(" << ZminG << "," << ZmaxG << ")";


// // nucSiteDistance is the min distance sites can be apart
dimensionedScalar lp("lp", dimensionSet(0, 1, 0, 0, 0),readScalar(PhaseFieldProperties.lookup("lp")));
dimensionedScalar lg("lg", dimensionSet(0, 1, 0, 0, 0),readScalar(PhaseFieldProperties.lookup("lg")));
scalar nucSiteDistance = nucDistFactor*max(lp.value(),lg.value());
Info << "Nucleation site distance: " << nucSiteDistance << endl;

//- RANDOM NUMBER GENERATORS
// Now I need to generate a List of N_max random points within this meshBox
// and then sync it across processors.
// Uniform distribution - for selection of nucleation sites
std::random_device rd;  // Will be used to obtain a seed for the random number engine
std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
std::uniform_real_distribution<> rand01(0.0, 1.0);


List<vector> nucCentres (N_max,vector(0,0,0));
List<bool> nucSet (N_max,false);
int nSet = 0;

scalar maxItersNucSet // how many interface thicknesses each nuc site should be from each other
(
    readScalar(nucleationDict.lookup("maxItersNucSet"))
);

Info << "Generating list of nucleation site locations" << endl;
// Do it in serial on the master processor, then sync to others
if (Pstream::master)
{
    forAll(nucCentres,i){
        bool nucOK = true;
        int iters = 0;

        vector nucC;

        do {

            // Assume nucOK is true
            nucOK = true;

            // Generate a random point within the domain
            scalar rX = rand01(gen);
            scalar rY = rand01(gen);
            scalar rZ = rand01(gen);

            scalar randX = XminG + rX*(XmaxG-XminG);
            scalar randY = YminG + rY*(YmaxG-YminG);
            scalar randZ = ZminG + rZ*(ZmaxG-ZminG);

            nucC = vector(randX,randY,randZ);
            nucCentres[i] = nucC;
            scalar R;
            
            // Need to check if this nucC conflicts with any others
            forAll(nucCentres,j){
                if (    (j != i) // dont compare to itself
                    &&  (nucSet[j] == true) // Just consider sites that have been already set
                    )
                {
                    // Get nucCentre that is being checked
                    vector nucNbr = nucCentres[j];

                    // get distance of existing nbr site from new nuc site
                    R = mag(nucC-nucNbr);
                    
                    if (R < nucSiteDistance) // too close to existing site
                    {
                        // Setting a site here will be invalid
                        nucOK = false;
                    }
                }
            }

            // This potential site has been checked against all other sites that have been set
            if (nucOK)
            {
                nucCentres[i] = nucC;
                nucSet[i] = true;
                nSet++;

                // Assign nucleation undercooling
                // forAll(meshC,cellI){
                //     scalar R = mag(nucC-meshC[cellI]);
                //     if (R < nucSiteDistance)
                //     {
                //         cellResidents[cellI] +=1;
                //     }
                // }
            }

            iters++;
            // Info << "i = " << i << "\titer: " << iters << "\tP = " <<  nucCentres[i] << endl;
            // int maxIters = 1e3;
            if (iters > maxItersNucSet)
            {
                FatalError << "Setting nucleation sites failed. " <<iters <<" iters taken (max = " << maxItersNucSet << ")" << abort(FatalError);
            }
        } while (nucOK == false); // keep trying until a valid site is found

        Info << "site[" << i << "] (" << 100.0*scalar(i)/scalar(nucCentres.size()) <<"%): " << nucC << "\titers: " << iters << endl;
    }
    Info << "Number of nucleations site set: " << nSet << " (N_max = " << N_max << ")" << endl;
    // Check if all nucleation sites have been set
    forAll(nucSet,k){
    	if (!nucSet[k])
    	{
    		FatalError << "Nuc site " << k << " not set!!" 
    		<<	"\tsetNucleationSite.H"
    		<< abort(FatalError);
    	}
    }
}

//- Send this list of nucCentres to each of the processors
Pstream::broadcast(nucCentres);
Pstream::broadcast(nucSet);


// Colour sites by processor number, just for interest
volScalarField procID
(
    IOobject
    (
        "procID",
        runTime.timeName(),
        mesh,
        // IOobject::READ_IF_PRESENT,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    Pstream::myProcNo(),
    zeroGradientFvPatchScalarField::typeName
);
procID.write();

volScalarField nucSiteProc
(
    IOobject
    (
        "nucSiteProc",
        runTime.timeName(),
        mesh,
        // IOobject::READ_IF_PRESENT,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    -1,
    zeroGradientFvPatchScalarField::typeName
);


int numNucSitesSet(0);
int numNucTuZero(0);
scalar totalTu(0);

//- List of nuc centres is now distributed to all processors
//	Need to find the cell centres that closest match these
if (nucCentres.size() > 0)
{
    Info << "nucCentres.size() = " << nucCentres.size() << endl;

	// Now I will assign each of the nucCentres using findCell
	List<label> nucCellIlist(nucCentres.size(),-1); // Initialise to -1
	forAll(nucCentres,i){
	    // find cell will return -1 if the cell isnt on the processor:
	    // get the nearest cell
	    label nucCellI = mesh.findCell(nucCentres[i]);
	    if (nucCellI != -1 && nucCellI < mesh.nCells()) // -1,>nCells if outside domain/processor
	    {
	        nucSite[nucCellI] = 1.0;
	        nucSiteProc[nucCellI] = Pstream::myProcNo();
	        nucCellIlist[i] = nucCellI;
	        Pout << "Site " << i << "\tfound by proc " << Pstream::myProcNo() << " at " << nucCentres[i] << endl;
	    }

	    // if (mag(nucCentres[i])<SMALL)
	    // {
	    //     FatalError << "Nucleation site unassigned" << abort(FatalError);
	    // }
	}

	// Pstream::scatter(dummy); // Sync
	nucSite.write();
	nucSiteProc.write();

	// Check for conflicts, just to make sure
	forAll(nucCentres,i){
	    vector parentNucC = nucCentres[i];

	    forAll(nucCentres,j){
	        if (i!=j) // Dont compare to itself
	        {
	            vector nbrNucC = nucCentres[j];
	            scalar R = mag(parentNucC - nbrNucC);
	            if (R < nucSiteDistance)
	            {
	                FatalError <<"Nucleation sites " << i << " and " << j << " are touching!!!!" << abort(FatalError);
	            }
	        }
	    }
	}
// } // Moving if to only set sites if N_max != 0

// int dummy;
// Pstream::scatter(dummy);
// nucSite.write();
// Info << "nucSite max: " << gMax(nucSite) << "\tmin: " << gMin(nucSite) << endl;
// abort();

Info << "Setting nucleation undercoolings for:\n"
<< "\tst dev:\t" << Tu_stdev.value() << endl
<< "\tmean:\t" << Tu_mean.value() << endl;

// int numNucSitesSet(0);
// int numNucTuZero(0);
// scalar totalTu(0);


// Gaussian distribution - for assignment of nucTu
std::default_random_engine generator;
std::normal_distribution<double> Tu_normal(Tu_mean.value(),Tu_stdev.value());

forAll(nucTu,i){
    if (nucSite[i] >= 0 )
    {
        scalar GaussTu = Tu_normal(generator);
        // scalar GaussTu = (GaussZ * Tu_stdev.value()) + Tu_mean.value();

        if (GaussTu < 0)
        {
            GaussTu = 0;
            numNucTuZero++;
        }

		Info << "nucSite[" << i << "]" << "\tTu: " << GaussTu << endl;


        // if (grainNum[i] < -0.5)
        // {
        //     nucSite[i] = 1.0; // 1 = active nuc site
        // }
        // else{
        //     nucSite[i] = 0.0; // 0 = inactive nuc site, will be activated once melted (xi < xiThreshold)
        // }
        
        nucTu[i] = GaussTu;

        totalTu += GaussTu;
        numNucSitesSet++;
    }
}
}

reduce(totalTu,sumOp<scalar>());
reduce(numNucSitesSet,sumOp<int>());

scalar avgNucTu(0);
if (numNucSitesSet > 0)
{
    avgNucTu = totalTu / numNucSitesSet;
}
Info << "numNucSitesSet: " << numNucSitesSet << "\tN_max: " << N_max << endl;
Info << "nucTu max: " << gMax(nucTu) << "\tmin: "<< gMin(nucTu) << "\tavg: " << avgNucTu<< "\tzero: " << numNucTuZero << endl;


Info << "\n\tSETTING OF NUCLEATION SITES COMPLETE" << nl << endl;


Switch stopToCheckNucSites
(
    nucleationDict.lookup("stopToCheckNucSites")
);
if (stopToCheckNucSites)
{
    nucTu.write();
    Info << "Stopping run to check nucleation sites (nucTu)" << endl;
    abort();
}
// nucTu.write();

// abort();
/*
Pstream::broadcast(dummy);

// Print range of x,y,z coordinates of all nuc sites
scalar minX, minY, minZ = GREAT;
scalar maxX, maxY, maxZ = SMALL;

forAll(nucCentres,i)
{
    minX = min(minX,nucCentres[i][0]);
    minY = min(minY,nucCentres[i][1]);
    minZ = min(minZ,nucCentres[i][2]);
    maxX = max(maxX,nucCentres[i][0]);
    maxY = max(maxY,nucCentres[i][1]);
    maxZ = max(maxZ,nucCentres[i][2]);
}

reduce(minX,minOp<scalar>());
reduce(minY,minOp<scalar>());
reduce(minZ,minOp<scalar>());
reduce(maxX,maxOp<scalar>());
reduce(maxY,maxOp<scalar>());
reduce(maxZ,maxOp<scalar>());

// Should print 

Info << "Nucleation site ranges:" << endl;
Info << "\tX: (" << minX << ", " << maxX << ")" << endl;
Info << "\tY: (" << minY << ", " << maxY << ")" << endl;
Info << "\tZ: (" << minZ << ", " << maxZ << ")" << endl;

Info << "\n\n\n\n\n\n\n" << endl;

// abort();

// Need to synchronise ni_fieldsNeeded
int nProcs = Pstream::nProcs();
int myProc = Pstream::myProcNo();

Info << "\nCarrying out reduce operation to sum up all number of ni_needed across "<< nProcs<<" processors" << endl;

int ni_fieldsNeeded_myProc = ni_fieldsNeeded;

Info << "before ni_fieldsNeeded: " << ni_fieldsNeeded << endl; 
reduce(ni_fieldsNeeded, sumOp<int>());
Info << "after ni_fieldsNeeded: " << ni_fieldsNeeded << endl;


for(int i = 0; i < nProcs; i++){
    if (i == myProc)
    {
        printf("Proc %d has %d sites and thinks the total is %d\n", myProc, ni_fieldsNeeded_myProc, ni_fieldsNeeded);
    }
}


Switch checkThatNucSitesExist
(
    nucleationDict.lookup("checkThatNucSitesExist")
);

if (checkThatNucSitesExist)
{
    if (gMax(nucTu) == 0)
    {
        FatalError << "No nucleation sites are set, check nucleation parameters and xiThreshold used in nucleationDict." << nl
        << "Domain could be already solid, if this is intentional, set checkThatNucSitesExist to false." << abort(FatalError);
    }
}

Info << "Number of potential nucleation sites: " << numNucSitesSet << endl;


// label numSitesToBeNuced(0);
numSitesToBeNuced = 0;

forAll(nucSite,i){
    if (nucSite[i] > -0.5){

        numSitesToBeNuced++;
    }
}


Info << "\n\nNo. potential nucleation sites: " << numSitesToBeNuced << nl << endl;

nucSite.write();
nucTu.write();

*/