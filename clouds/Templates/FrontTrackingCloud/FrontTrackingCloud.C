/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "FrontTrackingCloud.H"
//#include "IntegrationScheme.H"
#include "interpolation.H"
#include "subCycleTime.H"

//#include "InjectionModelList.H"
//#include "DispersionModel.H"
#include "FrontAnalyzingModel.H"
//#include "StochasticCollisionModel.H"
//#include "SurfaceFilmModel.H"
//#include "ModulationModel.H" //added ttw

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class CloudType>
void Foam::FrontTrackingCloud<CloudType>::setModels()
{
/*
    dispersionModel_.reset
    (
        DispersionModel<FrontTrackingCloud<CloudType> >::New
        (
            subModelProperties_,
            *this
        ).ptr()
    );
*/
    FrontAnalyzingModel_.reset
    (
        FrontAnalyzingModel<FrontTrackingCloud<CloudType> >::New
        (
            subModelProperties_,
            *this
        ).ptr()
    );
/*
    stochasticCollisionModel_.reset
    (
        StochasticCollisionModel<FrontTrackingCloud<CloudType> >::New
        (
            subModelProperties_,
            *this
        ).ptr()
    );

    surfaceFilmModel_.reset
    (
        SurfaceFilmModel<FrontTrackingCloud<CloudType> >::New
        (
            subModelProperties_,
            *this
        ).ptr()
    );

	modulationModel_.reset //added ttw
    (
        ModulationModel<FrontTrackingCloud<CloudType> >::New
        (
            subModelProperties_,
            *this
        ).ptr()
    );

    UIntegrator_.reset
    (
        vectorIntegrationScheme::New
        (
            "U",
            solution_.integrationSchemes()
        ).ptr()
    );
*/
}


template<class CloudType>
template<class TrackData>
void Foam::FrontTrackingCloud<CloudType>::solve(TrackData& td)
{
    if (solution_.steadyState())
    {
        td.cloud().storeState();

        td.cloud().preEvolve();

        evolveCloud(td);

        if (solution_.coupled())
        {
            td.cloud().relaxSources(td.cloud().cloudCopy());
        }
    }
    else
    {
        td.cloud().preEvolve();

        evolveCloud(td);

        if (solution_.coupled())
        {
            td.cloud().scaleSources();
        }
    }

    td.cloud().info();

    //Std.cloud().postEvolve();

    if (solution_.steadyState())
    {
        td.cloud().restoreState();
    }
}


template<class CloudType>
void Foam::FrontTrackingCloud<CloudType>::buildCellOccupancy()
{
    if (cellOccupancyPtr_.empty())
    {
        cellOccupancyPtr_.reset
        (
            new List<DynamicList<parcelType*> >(mesh_.nCells())
        );
    }
    else if (cellOccupancyPtr_().size() != mesh_.nCells())
    {
        // If the size of the mesh has changed, reset the
        // cellOccupancy size

        cellOccupancyPtr_().setSize(mesh_.nCells());
    }

    List<DynamicList<parcelType*> >& cellOccupancy = cellOccupancyPtr_();

    forAll(cellOccupancy, cO)
    {
        cellOccupancy[cO].clear();
    }

    forAllIter(typename FrontTrackingCloud<CloudType>, *this, iter)
    {
        cellOccupancy[iter().cell()].append(&iter());
    }
}


template<class CloudType>
void Foam::FrontTrackingCloud<CloudType>::updateCellOccupancy()
{
    // Only build the cellOccupancy if the pointer is set, i.e. it has
    // been requested before.

    if (cellOccupancyPtr_.valid())
    {
        buildCellOccupancy();
    }
}


template<class CloudType>
template<class TrackData>
void Foam::FrontTrackingCloud<CloudType>::evolveCloud(TrackData& td)
{
    if (solution_.coupled())
    {
                // ***you must use foam error to note choosing the coupled option. 
        td.cloud().resetSourceTerms();
    }

    if (solution_.transient())
    {
                // ***you must use foam error to note choosing the transient option. 
        
        //label preInjectionSize = this->size();

        //this->surfaceFilm().inject(td);

        // Update the cellOccupancy if the size of the cloud has changed
        // during the injection.
        /*
        if (preInjectionSize != this->size())
        {
            updateCellOccupancy();
            preInjectionSize = this->size();
        }
        */

        //injectors_.inject(td);

		
        // Assume that motion will update the cellOccupancy as necessary
        // before it is required.

        FrontAnalyzing().update(); 

        td.cloud().motion(td);

            Info << ' ' << "I am after motion call" << endl;

                // ***maybe you can use it.
		//td.cloud().preCollisionUpdates(td); //added1

        //stochasticCollision().update(solution_.trackTime()); 

        // it is seem that it is better to pass td in function.

    }
/*
    else
    {
//        this->surfaceFilm().injectSteadyState(td);

        injectors_.injectSteadyState(td, solution_.trackTime());

        td.part() = TrackData::tpLinearTrack;
        CloudType::move(td,  solution_.trackTime());
    }
*/

}


template<class CloudType>
void Foam::FrontTrackingCloud<CloudType>::postEvolve()
{
    Info<< endl;

    if (debug)
    {
        this->writePositions();
    }

    //this->dispersion().cacheFields(false);

    //forces_.cacheFields(false);

    //functions_.postEvolve();

    solution_.nextIter();

    if (this->db().time().outputTime())
    {
        outputProperties_.writeObject
        (
            IOstream::ASCII,
            IOstream::currentVersion,
            this->db().time().writeCompression()
        );
    }
}


template<class CloudType>
void Foam::FrontTrackingCloud<CloudType>::cloudReset(FrontTrackingCloud<CloudType>& c)
{
    CloudType::cloudReset(c);

    rndGen_ = c.rndGen_;

    //forces_.transfer(c.forces_);

    //functions_.transfer(c.functions_);

    //injectors_.transfer(c.injectors_);

    //dispersionModel_.reset(c.dispersionModel_.ptr());
    FrontAnalyzingModel_.reset(c.FrontAnalyzingModel_.ptr());
    //stochasticCollisionModel_.reset(c.stochasticCollisionModel_.ptr());
    //surfaceFilmModel_.reset(c.surfaceFilmModel_.ptr());
    //modulationModel_.reset(c.modulationModel_.ptr()); //added ttw

    //UIntegrator_.reset(c.UIntegrator_.ptr());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::FrontTrackingCloud<CloudType>::FrontTrackingCloud
(
    const word& cloudName,
    const volScalarField& rho,
    const volVectorField& U,
    const volScalarField& mu,
    const volScalarField& p,//Jafari added
    const dimensionedVector& g,
    bool readFields
)
:
    CloudType(rho.mesh(), cloudName, false),
    frontTrackingCloud(),
    cloudCopyPtr_(NULL),
    mesh_(rho.mesh()),
    particleProperties_
    (
        IOobject
        (
            cloudName + "Properties",
            rho.mesh().time().constant(),
            rho.mesh(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    outputProperties_
    (
        IOobject
        (
            cloudName + "OutputProperties",
            mesh_.time().timeName(),
            "uniform"/cloud::prefix/cloudName,
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    ),
    solution_(mesh_, particleProperties_.subDict("solution")),
    constProps_(particleProperties_),
    subModelProperties_
    (
        particleProperties_.subOrEmptyDict("subModels", solution_.active())
    ),
    rndGen_
    (
        label(0),
        solution_.steadyState() ?
        particleProperties_.lookupOrDefault<label>("randomSampleSize", 100000)
      : -1
    ),
    cellOccupancyPtr_(),
    cellLengthScale_(cbrt(mesh_.V())),
    rho_(rho),
    U_(U),
    mu_(mu),
    p_(p),//Jafari added
    g_(g),
    pAmbient_(0.0),
/*
    forces_
    (
        *this,
        mesh_,
        subModelProperties_.subOrEmptyDict
        (
            "particleForces",
            solution_.active()
        ),
        solution_.active()
    ),
*/
/*
    functions_
    (
        *this,
        particleProperties_.subOrEmptyDict("cloudFunctions"),
        solution_.active()
    ),
*/
/*
    injectors_
    (
        subModelProperties_.subOrEmptyDict("injectionModels"),
        *this
    ),
    dispersionModel_(NULL),
*/
    FrontAnalyzingModel_(NULL),
/*
    stochasticCollisionModel_(NULL),
    surfaceFilmModel_(NULL),
*/
/*
	modulationModel_(NULL), //added ttw
	preMotStat1Comp_(false), //added1
	preMotStat2Comp_(false),
	preColStat1Comp_(false),
	preColStat2Comp_(false),
*/
/*
    UIntegrator_(NULL),
*/

/*
    UTrans_
    (
        new DimensionedField<vector, volMesh>
        (
            IOobject
            (
                this->name() + ":UTrans",
                this->db().time().timeName(),
                this->db(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedVector("zero", dimMass*dimVelocity, vector::zero)
        )
    ),
    UCoeff_
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                this->name() + ":UCoeff",
                this->db().time().timeName(),
                this->db(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero",  dimMass, 0.0)
        )
    ),
    RTrans_ //added ttw
    (
        new DimensionedField<symmTensor, volMesh>
        (
            IOobject
            (
                this->name() + ":RTrans",
                this->db().time().timeName(),
                this->db(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedSymmTensor("zero", dimMass*dimVelocity*dimVelocity, symmTensor::zero)
        )
    )
*/
    densityFromFT_
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                this->name() + ":densityFromFT",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", dimDensity, 0.0)
        )
    ),
    viscosityFromFT_
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                this->name() + ":viscosityFromFT",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", dimMass/dimLength/dimTime, 0.0)
        )
    ),
    sTensionForceFromFT_
    (
        new DimensionedField<vector, volMesh>
        (
            IOobject
            (
                this->name() + ":sTensionForceFromFT",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedVector("sTensionForceFromFT", dimForce/dimVolume,  vector::zero)
        )
    ),
    dPressureFromFT_
    (
        new DimensionedField<vector, volMesh>
        (
            IOobject
            (
                this->name() + ":dPressureFromFT",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedVector("dPressureFromFT", dimForce/dimVolume,  vector::zero)
        )
    ),
    pressureJumpAtTheInterfaceFromFT_
    (
        new DimensionedField<vector, volMesh>
        (
            IOobject
            (
                this->name() + ":pressureJumpAtTheInterfaceFromFT",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedVector("pressureJumpAtTheInterfaceFromFT", dimForce/dimVolume,  vector::zero)
        )
    )
{
    if (solution_.active())
    {
        setModels();

        if (readFields)
        {
            parcelType::readFields(*this);
        }
    }

    if (solution_.resetSourcesOnStartup())
    {
        // **** maybe this can be used
        resetSourceTerms();
    }
}


template<class CloudType>
Foam::FrontTrackingCloud<CloudType>::FrontTrackingCloud
(
    FrontTrackingCloud<CloudType>& c,
    const word& name
)
:
    CloudType(c.mesh_, name, c),
    frontTrackingCloud(),
    cloudCopyPtr_(NULL),
    mesh_(c.mesh_),
    particleProperties_(c.particleProperties_),
    outputProperties_(c.outputProperties_),
    solution_(c.solution_),
    constProps_(c.constProps_),
    subModelProperties_(c.subModelProperties_),
    rndGen_(c.rndGen_, true),
    cellOccupancyPtr_(NULL),
    cellLengthScale_(c.cellLengthScale_),
    rho_(c.rho_),
    U_(c.U_),
    mu_(c.mu_),
    p_(c.p_),//Jafari added
    g_(c.g_),
    pAmbient_(c.pAmbient_),
    //forces_(c.forces_),
    //functions_(c.functions_),
    //injectors_(c.injectors_),
    //dispersionModel_(c.dispersionModel_->clone()),
    FrontAnalyzingModel_(c.FrontAnalyzingModel_->clone()),
    //stochasticCollisionModel_(c.stochasticCollisionModel_->clone()),
    //surfaceFilmModel_(c.surfaceFilmModel_->clone()),
	//modulationModel_(c.modulationModel_->clone()), //added ttw
	//preMotStat1Comp_(c.preMotStat1Comp_), //added1
	//preMotStat2Comp_(c.preMotStat2Comp_), //added1
	//preColStat1Comp_(c.preColStat1Comp_),
	//preColStat2Comp_(c.preColStat2Comp_),
    //UIntegrator_(c.UIntegrator_->clone()),
/*
    UTrans_
    (
        new DimensionedField<vector, volMesh>
        (
            IOobject
            (
                this->name() + ":UTrans",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            c.UTrans_()
        )
    ),
    UCoeff_
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                name + ":UCoeff",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            c.UCoeff_()
        )
    ),
    RTrans_ //added ttw
    (
        new DimensionedField<symmTensor, volMesh>
        (
            IOobject
            (
                this->name() + ":RTrans",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            c.RTrans_()
        )
    )
*/
    densityFromFT_
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                name + ":densityFromFT",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            c.densityFromFT_()
        )
    ),
    viscosityFromFT_
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                name + ":viscosityFromFT",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            c.viscosityFromFT_()
        )
    ),
    sTensionForceFromFT_
    (
        new DimensionedField<vector, volMesh>
        (
            IOobject
            (
                this->name() + ":sTensionForceFromFT",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            c.sTensionForceFromFT_()
        )
    ),
    dPressureFromFT_
    (
        new DimensionedField<vector, volMesh>
        (
            IOobject
            (
                this->name() + ":dPressureFromFT",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            c.dPressureFromFT_()
        )
    ),
    pressureJumpAtTheInterfaceFromFT_
    (
        new DimensionedField<vector, volMesh>
        (
            IOobject
            (
                this->name() + ":pressureJumpAtTheInterfaceFromFT",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            c.pressureJumpAtTheInterfaceFromFT_()
        )
    )
{}


template<class CloudType>
Foam::FrontTrackingCloud<CloudType>::FrontTrackingCloud
(
    const fvMesh& mesh,
    const word& name,
    const FrontTrackingCloud<CloudType>& c
)
:
    CloudType(mesh, name, IDLList<parcelType>()),
    frontTrackingCloud(),
    cloudCopyPtr_(NULL),
    mesh_(mesh),
    particleProperties_
    (
        IOobject
        (
            name + "Properties",
            mesh.time().constant(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        )
    ),
    outputProperties_
    (
        IOobject
        (
            name + "OutputProperties",
            mesh_.time().timeName(),
            "uniform"/cloud::prefix/name,
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        )
    ),
    solution_(mesh),
    constProps_(),
    subModelProperties_(dictionary::null),
    rndGen_(0, 0),
    cellOccupancyPtr_(NULL),
    cellLengthScale_(c.cellLengthScale_),
    rho_(c.rho_),
    U_(c.U_),
    mu_(c.mu_),
    p_(c.p_),//Jafari added
    g_(c.g_),
    pAmbient_(c.pAmbient_),
    //forces_(*this, mesh),
    //functions_(*this),
    //injectors_(*this),
    //dispersionModel_(NULL),
    FrontAnalyzingModel_(NULL),
    //stochasticCollisionModel_(NULL),
    //surfaceFilmModel_(NULL),
	//modulationModel_(NULL), //added ttw
	//preMotStat1Comp_(false), //added1
	//preMotStat2Comp_(false),
	//preColStat1Comp_(false), //added1
	//preColStat2Comp_(false),
    //UIntegrator_(NULL),
    //UTrans_(NULL),
    //UCoeff_(NULL),
	//RTrans_(NULL) //added ttw
    densityFromFT_(NULL),
    viscosityFromFT_(NULL),
    sTensionForceFromFT_(NULL),
    dPressureFromFT_(NULL),
    pressureJumpAtTheInterfaceFromFT_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::FrontTrackingCloud<CloudType>::~FrontTrackingCloud()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
bool Foam::FrontTrackingCloud<CloudType>::hasWallImpactDistance() const
{
    return true;
}


template<class CloudType>
void Foam::FrontTrackingCloud<CloudType>::setParcelThermoProperties
(
    parcelType& parcel,
    const scalar lagrangianDt
)
{
    parcel.rho() = constProps_.rho0();
}

// *** maybe is can not be benefit
template<class CloudType>
void Foam::FrontTrackingCloud<CloudType>::checkParcelProperties
(
    parcelType& parcel,
    const scalar lagrangianDt,
    const bool fullyDescribed
)
{
    const scalar carrierDt = mesh_.time().deltaTValue();
    parcel.stepFraction() = (carrierDt - lagrangianDt)/carrierDt;

    if (parcel.typeId() == -1)
    {
        parcel.typeId() = constProps_.parcelTypeId();
    }
}


template<class CloudType>
void Foam::FrontTrackingCloud<CloudType>::storeState()
{
    cloudCopyPtr_.reset
    (
        static_cast<FrontTrackingCloud<CloudType>*>
        (
            clone(this->name() + "Copy").ptr()
        )
    );
}


template<class CloudType>
void Foam::FrontTrackingCloud<CloudType>::restoreState()
{
    cloudReset(cloudCopyPtr_());
    cloudCopyPtr_.clear();
}

//*** maybe it can be used again
template<class CloudType>
void Foam::FrontTrackingCloud<CloudType>::resetSourceTerms()
{
    //UTrans().field() = vector::zero;
    //UCoeff().field() = 0.0;
	//RTrans().field() = symmTensor::zero; //added ttw
}

//*** maybe it can be used again
template<class CloudType>
template<class Type>
void Foam::FrontTrackingCloud<CloudType>::relax
(
    DimensionedField<Type, volMesh>& field,
    const DimensionedField<Type, volMesh>& field0,
    const word& name
) const
{
    //const scalar coeff = solution_.relaxCoeff(name);
    //field = field0 + coeff*(field - field0);
}

//*** maybe it can be used again
template<class CloudType>
template<class Type>
void Foam::FrontTrackingCloud<CloudType>::scale
(
    DimensionedField<Type, volMesh>& field,
    const word& name
) const
{
   // const scalar coeff = solution_.relaxCoeff(name);
    //field *= coeff;
}

//*** maybe it can be used again
template<class CloudType>
void Foam::FrontTrackingCloud<CloudType>::relaxSources
(
    const FrontTrackingCloud<CloudType>& cloudOldTime
)
{
   // this->relax(UTrans_(), cloudOldTime.UTrans(), "U");
   // this->relax(UCoeff_(), cloudOldTime.UCoeff(), "U");
	//this->relax(RTrans_(), cloudOldTime.RTrans(), "R"); //added ttw
}

//*** maybe it can be used again
template<class CloudType>
void Foam::FrontTrackingCloud<CloudType>::scaleSources()
{
    //this->scale(UTrans_(), "U");
    //this->scale(UCoeff_(), "U");
	//this->scale(RTrans_(), "R"); //added ttw
}


//*** maybe it can be used again
template<class CloudType>
void Foam::FrontTrackingCloud<CloudType>::preEvolve()
{
    // force calculaion of mesh dimensions - needed for parallel runs
    // with topology change due to lazy evaluation of valid mesh dimensions
    label nGeometricD = mesh_.nGeometricD();

    Info<< "\nSolving " << nGeometricD << "-D cloud " << this->name() << endl;

    //this->dispersion().cacheFields(true);
    //forces_.cacheFields(true);
    updateCellOccupancy();

    Info<< "\nAFter...........................Solving " << endl;
   pAmbient_ = constProps_.dict().template
        lookupOrDefault<scalar>("pAmbient", pAmbient_);

    //functions_.preEvolve();
}


template<class CloudType>
void Foam::FrontTrackingCloud<CloudType>::evolve()
{
    if (solution_.canEvolve())
   {
        typename parcelType::template
            TrackingData<FrontTrackingCloud<CloudType> > td(*this);

		//td.updateParticleStatistics(*this); //added1

        solve(td);
    Info<< "Cloud: " << this->name() << endl;
   }
}


template<class CloudType>
template<class TrackData>
void Foam::FrontTrackingCloud<CloudType>::motion(TrackData& td)
{
    td.part() = TrackData::tpLinearTrack;
    CloudType::move(td, solution_.trackTime());
            Info << ' ' << "I am after motion/move call" << endl;

    updateCellOccupancy();
}


template<class CloudType>
template<class TrackData>
void Foam::FrontTrackingCloud<CloudType>::preCollisionUpdates(TrackData& td)
{
	td.updates() = TrackData::statCollision; //added1
	td.updateParticleStatistics(); //(*this);
}


template<class CloudType>
void Foam::FrontTrackingCloud<CloudType>::patchData
(
    const parcelType& p,
    const polyPatch& pp,
    const scalar trackFraction,
    const tetIndices& tetIs,
    vector& nw,
    vector& Up
) const
{
    label patchI = pp.index();
    label patchFaceI = pp.whichFace(p.face());

    vector n = tetIs.faceTri(mesh_).normal();
    n /= mag(n);

    vector U = U_.boundaryField()[patchI][patchFaceI];

    // Unless the face is rotating, the required normal is n;
    nw = n;

    if (!mesh_.moving())
    {
        // Only wall patches may have a non-zero wall velocity from
        // the velocity field when the mesh is not moving.

        if (isA<wallPolyPatch>(pp))
        {
            Up = U;
        }
        else
        {
            Up = vector::zero;
        }
    }
    else
    {
        vector U00 = U_.oldTime().boundaryField()[patchI][patchFaceI];

        vector n00 = tetIs.oldFaceTri(mesh_).normal();

        // Difference in normal over timestep
        vector dn = vector::zero;

        if (mag(n00) > SMALL)
        {
            // If the old normal is zero (for example in layer
            // addition) then use the current normal, meaning that the
            // motion can only be translational, and dn remains zero,
            // otherwise, calculate dn:

            n00 /= mag(n00);

            dn = n - n00;
        }

        // Total fraction through the timestep of the motion,
        // including stepFraction before the current tracking step
        // and the current trackFraction
        // i.e.
        // let s = stepFraction, t = trackFraction
        // Motion of x in time:
        // |-----------------|---------|---------|
        // x00               x0        xi        x
        //
        // where xi is the correct value of x at the required
        // tracking instant.
        //
        // x0 = x00 + s*(x - x00) = s*x + (1 - s)*x00
        //
        // i.e. the motion covered by previous tracking portions
        // within this timestep, and
        //
        // xi = x0 + t*(x - x0)
        //    = t*x + (1 - t)*x0
        //    = t*x + (1 - t)*(s*x + (1 - s)*x00)
        //    = (s + t - s*t)*x + (1 - (s + t - s*t))*x00
        //
        // let m = (s + t - s*t)
        //
        // xi = m*x + (1 - m)*x00 = x00 + m*(x - x00);
        //
        // In the same form as before.

        scalar m =
            p.stepFraction()
          + trackFraction
          - (p.stepFraction()*trackFraction);

        // When the mesh is moving, the velocity field on wall patches
        // will contain the velocity associated with the motion of the
        // mesh, in which case it is interpolated in time using m.
        // For other patches the face velocity will need to be
        // reconstructed from the face centre motion.

        const vector& Cf = mesh_.faceCentres()[p.face()];

        vector Cf00 = mesh_.faces()[p.face()].centre(mesh_.oldPoints());

        if (isA<wallPolyPatch>(pp))
        {
            Up = U00 + m*(U - U00);
        }
        else
        {
            Up = (Cf - Cf00)/mesh_.time().deltaTValue();
        }

        if (mag(dn) > SMALL)
        {
            // Rotational motion, nw requires interpolation and a
            // rotational velocity around face centre correction to Up
            // is required.

            nw = n00 + m*dn;

            // Cf at tracking instant
            vector Cfi = Cf00 + m*(Cf - Cf00);

            // Normal vector cross product
            vector omega = (n00 ^ n);

            scalar magOmega = mag(omega);

            // magOmega = sin(angle between unit normals)
            // Normalise omega vector by magOmega, then multiply by
            // angle/dt to give the correct angular velocity vector.
            omega *= Foam::asin(magOmega)/(magOmega*mesh_.time().deltaTValue());

            // Project position onto face and calculate this position
            // relative to the face centre.
            vector facePos =
                p.position()
              - ((p.position() - Cfi) & nw)*nw
              - Cfi;

            Up += (omega ^ facePos);
        }

        // No further action is required if the motion is
        // translational only, nw and Up have already been set.
    }
}


template<class CloudType>
void Foam::FrontTrackingCloud<CloudType>::updateMesh()
{
    updateCellOccupancy();
    //injectors_.updateMesh();
    cellLengthScale_ = cbrt(mesh_.V());
}


template<class CloudType>
void Foam::FrontTrackingCloud<CloudType>::autoMap(const mapPolyMesh& mapper)
{
    typedef typename particle::TrackingData<FrontTrackingCloud<CloudType> > tdType;

    tdType td(*this);

    Cloud<parcelType>::template autoMap<tdType>(td, mapper);

    updateMesh();
}


template<class CloudType>
void Foam::FrontTrackingCloud<CloudType>::info()
{
    //vector linearMomentum = linearMomentumOfSystem();
    //reduce(linearMomentum, sumOp<vector>());

    //scalar linearKineticEnergy = linearKineticEnergyOfSystem();
    //reduce(linearKineticEnergy, sumOp<scalar>());

/*
    Info<< "Cloud: " << this->name() << nl
        << "    Current number of parcels       = "
        << returnReduce(this->size(), sumOp<label>()) << nl
        << "    Current mass in system          = "
        << returnReduce(massInSystem(), sumOp<scalar>()) << nl
        << "    Linear momentum                 = "
        << linearMomentum << nl
        << "   |Linear momentum|                = "
        << mag(linearMomentum) << nl
        << "    Linear kinetic energy           = "
        << linearKineticEnergy << nl;
*/

    //injectors_.info(Info);
    //this->surfaceFilm().info(Info);
    //this->FrontAnalyzing().info(Info);
}


// ************************************************************************* //
