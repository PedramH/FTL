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

#include "KinematicRotCloud.H"
#include "KinematicRotParcel.H"
//#include "fvcCurl.H" //OmegaC memory

//#include "DeterminesticCollisionModel.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

//OmegaC memory
/*template<class CloudType>
void Foam::KinematicRotCloud<CloudType>::cacheFields(const bool store)
{
    static word fName("curlUcDt");

    bool fieldExists = this->mesh().template foundObject<volVectorField>(fName);

    if (store)
    {
        if (!fieldExists)
        {
            const volVectorField& Uc = this->mesh().template
                lookupObject<volVectorField>(UName_);

            curlUPtr_ = new volVectorField(fName, fvc::curl(Uc));

            curlUPtr->store();
        }
    }
    else
    {
        if (fieldExists)
        {
            const volVectorField& curlUc = this->mesh().template
                lookupObject<volVectorField>(fName);

            const_cast<volVectorField&>(curlUc).checkOut();
        }
    }
}*/


template<class CloudType>
void Foam::KinematicRotCloud<CloudType>::setModels()
{
    /*determinesticCollisionModel_.reset
    (
        DeterminesticCollisionModel<KinematicRotCloud<CloudType> >::New
        (
            this->subModelProperties(),
            *this
        ).ptr()
    );*/

    OmegaIntegrator_.reset
    (
        vectorIntegrationScheme::New
        (
            "Omega",
            this->solution().integrationSchemes()
        ).ptr()
    );
}


template<class CloudType>
void Foam::KinematicRotCloud<CloudType>::cloudReset(KinematicRotCloud<CloudType>& c)
{
    CloudType::cloudReset(c);

    torques_.transfer(c.torques_); 
	//determinesticCollisionModel_.reset(c.determinesticCollisionModel_.ptr());

    OmegaIntegrator_.reset(c.OmegaIntegrator_.ptr());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::KinematicRotCloud<CloudType>::KinematicRotCloud
(
    const word& cloudName,
    const volScalarField& rho,
    const volVectorField& U,
    const volScalarField& mu,
    const dimensionedVector& g,
    bool readFields
)
:
	CloudType(cloudName, rho, U, mu, g, false),
    kinematicRotCloud(),
    cloudCopyPtr_(NULL),
    constProps_(this->particleProperties()),
    torques_ 
    (
        *this,
        this->mesh(),
        this->subModelProperties().subOrEmptyDict
        (
            "particleTorques",
            this->solution().active()
        ),
        this->solution().active()
    ),
	//curlUPtr_(NULL), //OmegaC memory
    //determinesticCollisionModel_(NULL),
	//stat2Comp_(false),
    OmegaIntegrator_(NULL),
	OmegaTrans_
    (
        new DimensionedField<vector, volMesh>
        (
            IOobject
            (
                this->name() + ":OmegaTrans",
                this->db().time().timeName(),
                this->db(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            this->mesh(),
            dimensionedVector("zero", dimForce*dimLength*dimTime, vector::zero)
        )
    )
    /*OmegaCoeff_
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                this->name() + ":OmegaCoeff",
                this->db().time().timeName(),
                this->db(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            this->mesh(),
            dimensionedScalar("zero",  dimMass*dimLength*dimLength, 0.0)
        )
    )*/ 
{
    if (this->solution().active())
    {
        setModels();

        if (readFields)
        {
            parcelType::readFields(*this);
        }
    }

    if (this->solution().resetSourcesOnStartup())
    {
        resetSourceTerms();
    }
}


template<class CloudType>
Foam::KinematicRotCloud<CloudType>::KinematicRotCloud
(
    KinematicRotCloud<CloudType>& c,
    const word& name
)
:
    CloudType(c, name),
    kinematicRotCloud(),
    cloudCopyPtr_(NULL),
	constProps_(c.constProps_),
    torques_(c.torques_),
	//curlUPtr_(c.curlUPtr_), //OmegaC memory
	//determinesticCollisionModel_(c.determinesticCollisionModel_->clone()),
	//stat2Comp_(c.stat2Comp_),
    OmegaIntegrator_(c.OmegaIntegrator_->clone()),
	OmegaTrans_
    (
        new DimensionedField<vector, volMesh>
        (
            IOobject
            (
                this->name() + ":OmegaTrans",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            c.OmegaTrans_()
        )
    )
    /*OmegaCoeff_
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                name + ":OmegaCoeff",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            c.OmegaCoeff_()
        )
    )*/
{}


template<class CloudType>
Foam::KinematicRotCloud<CloudType>::KinematicRotCloud
(
    const fvMesh& mesh,
    const word& name,
    const KinematicRotCloud<CloudType>& c
)
:
    CloudType(mesh, name, c),
    kinematicRotCloud(),
    cloudCopyPtr_(NULL),
	constProps_(),
    torques_(*this, mesh),
	//curlUPtr_(NULL), //OmegaC memory
    //determinesticCollisionModel_(NULL),
    //OmegaCoeff_(NULL), 
	//stat2Comp_(false),
    OmegaIntegrator_(NULL),
	OmegaTrans_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::KinematicRotCloud<CloudType>::~KinematicRotCloud()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::KinematicRotCloud<CloudType>::setParcelThermoProperties
(
    parcelType& parcel,
    const scalar lagrangianDt
)
{
    CloudType::setParcelThermoProperties(parcel, lagrangianDt);

    parcel.Omega() = constProps_.Omega0(); //for injection and surface film models
}


template<class CloudType>
void Foam::KinematicRotCloud<CloudType>::checkParcelProperties
(
    parcelType& parcel,
    const scalar lagrangianDt,
    const bool fullyDescribed
)
{
    CloudType::checkParcelProperties(parcel, lagrangianDt, fullyDescribed);
}


template<class CloudType>
void Foam::KinematicRotCloud<CloudType>::storeState()
{
    cloudCopyPtr_.reset
    (
        static_cast<KinematicRotCloud<CloudType>*>
        (
            clone(this->name() + "Copy").ptr()
        )
    );
}


template<class CloudType>
void Foam::KinematicRotCloud<CloudType>::restoreState()
{
    cloudReset(cloudCopyPtr_());
    cloudCopyPtr_.clear();
}


template<class CloudType>
void Foam::KinematicRotCloud<CloudType>::resetSourceTerms()
{
    CloudType::resetSourceTerms();
    OmegaTrans_->field() = vector::zero;
    //OmegaCoeff_->field() = 0.0;
}


template<class CloudType>
void Foam::KinematicRotCloud<CloudType>::relaxSources
(
    const KinematicRotCloud<CloudType>& cloudOldTime
)
{
    CloudType::relaxSources(cloudOldTime);

	this->relax(OmegaTrans_(), cloudOldTime.UTrans(), "Omega");
    //this->relax(OmegaCoeff_(), cloudOldTime.UCoeff(), "Omega");
}


template<class CloudType>
void Foam::KinematicRotCloud<CloudType>::scaleSources()
{
    CloudType::scaleSources();

    this->scale(OmegaTrans_(), "Omega");
    //this->scale(OmegaCoeff_(), "Omega");
}


template<class CloudType>
void Foam::KinematicRotCloud<CloudType>::preEvolve()
{
    CloudType::preEvolve();

    torques_.cacheFields(true); 

	//cacheFields(true); //OmegaC memory
}


template<class CloudType>
void Foam::KinematicRotCloud<CloudType>::evolve()
{
    if (this->solution().canEvolve())
    {
        typename parcelType::template
            TrackingData<KinematicRotCloud<CloudType> > td(*this);	
	
		//td.updateParticleStatistics(*this);    //added1

        this->solve(td);
    }
}


template<class CloudType>
template<class TrackData>
void Foam::KinematicRotCloud<CloudType>::motion(TrackData& td)
{
    CloudType::motion(td);
}


/*template<class CloudType> //no difference in commenting this function //added1
template<class TrackData>
void Foam::KinematicRotCloud<CloudType>::preCollisionUpdates(TrackData& td)
{
	CloudType::preCollisionUpdates(td);
}*/


template<class CloudType>
void Foam::KinematicRotCloud<CloudType>::postEvolve()
{
    CloudType::postEvolve();

    torques_.cacheFields(false);

	//cacheFields(false); //OmegaC memory
}


template<class CloudType>
void Foam::KinematicRotCloud<CloudType>::autoMap(const mapPolyMesh& mapper)
{
    typedef typename particle::TrackingData<KinematicRotCloud<CloudType> > tdType;

    tdType td(*this);

    Cloud<parcelType>::template autoMap<tdType>(td, mapper);

    this->updateMesh();
}


template<class CloudType>
void Foam::KinematicRotCloud<CloudType>::info()
{
    CloudType::info();

    Info<< "    Angualr velocity min/max             = " << Omegamin() << ", " << Omegamax()
        << endl;
}


// ************************************************************************* //
