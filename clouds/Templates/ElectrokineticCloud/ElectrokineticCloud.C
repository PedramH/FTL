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

#include "ElectrokineticCloud.H"
#include "ElectrokineticParcel.H"
//#include "fvcCurl.H" //OmegaC memory

//#include "DeterminesticCollisionModel.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class CloudType>
void Foam::ElectrokineticCloud<CloudType>::setModels()
{}


template<class CloudType>
void Foam::ElectrokineticCloud<CloudType>::cloudReset(ElectrokineticCloud<CloudType>& c)
{
    CloudType::cloudReset(c);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ElectrokineticCloud<CloudType>::ElectrokineticCloud
(
    const word& cloudName,
    const volScalarField& rho,
    const volVectorField& U,
    const volScalarField& mu,
	const volVectorField& E,
    const dimensionedVector& g,
    bool readFields
)
:
	CloudType(cloudName, rho, U, mu, g, false),
    electrokineticCloud(),
    cloudCopyPtr_(NULL),
    constProps_(this->particleProperties()),
    E_(E)
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
Foam::ElectrokineticCloud<CloudType>::ElectrokineticCloud
(
    ElectrokineticCloud<CloudType>& c,
    const word& name
)
:
    CloudType(c, name),
    electrokineticCloud(),
    cloudCopyPtr_(NULL),
	constProps_(c.constProps_),
    E_(c.E_)
{}


template<class CloudType>
Foam::ElectrokineticCloud<CloudType>::ElectrokineticCloud
(
    const fvMesh& mesh,
    const word& name,
    const ElectrokineticCloud<CloudType>& c
)
:
    CloudType(mesh, name, c),
    electrokineticCloud(),
    cloudCopyPtr_(NULL),
	constProps_(),
    E_(c.E_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ElectrokineticCloud<CloudType>::~ElectrokineticCloud()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::ElectrokineticCloud<CloudType>::setParcelThermoProperties
(
    parcelType& parcel,
    const scalar lagrangianDt
)
{
    CloudType::setParcelThermoProperties(parcel, lagrangianDt);

    parcel.z() = constProps_.z0(); //for injection and surface film models
}


template<class CloudType>
void Foam::ElectrokineticCloud<CloudType>::checkParcelProperties
(
    parcelType& parcel,
    const scalar lagrangianDt,
    const bool fullyDescribed
)
{
    CloudType::checkParcelProperties(parcel, lagrangianDt, fullyDescribed);
}


template<class CloudType>
void Foam::ElectrokineticCloud<CloudType>::storeState()
{
    cloudCopyPtr_.reset
    (
        static_cast<ElectrokineticCloud<CloudType>*>
        (
            clone(this->name() + "Copy").ptr()
        )
    );
}


template<class CloudType>
void Foam::ElectrokineticCloud<CloudType>::restoreState()
{
    cloudReset(cloudCopyPtr_());
    cloudCopyPtr_.clear();
}


template<class CloudType>
void Foam::ElectrokineticCloud<CloudType>::resetSourceTerms()
{
    CloudType::resetSourceTerms();
}


template<class CloudType>
void Foam::ElectrokineticCloud<CloudType>::relaxSources
(
    const ElectrokineticCloud<CloudType>& cloudOldTime
)
{
    CloudType::relaxSources(cloudOldTime);
}


template<class CloudType>
void Foam::ElectrokineticCloud<CloudType>::scaleSources()
{
    CloudType::scaleSources();
}


template<class CloudType>
void Foam::ElectrokineticCloud<CloudType>::preEvolve()
{
    CloudType::preEvolve();
}


template<class CloudType>
void Foam::ElectrokineticCloud<CloudType>::evolve()
{
    if (this->solution().canEvolve())
    {
        typename parcelType::template
            TrackingData<ElectrokineticCloud<CloudType> > td(*this);	
	
		//td.updateParticleStatistics(*this); //added1   

        this->solve(td);
    }
}


template<class CloudType>
template<class TrackData>
void Foam::ElectrokineticCloud<CloudType>::motion(TrackData& td)
{
    CloudType::motion(td);
}


/*template<class CloudType>
template<class TrackData>
void Foam::ElectrokineticCloud<CloudType>::preCollisionUpdates(TrackData& td)
{
	CloudType::preCollisionUpdates(td); //added1
}*/


template<class CloudType>
void Foam::ElectrokineticCloud<CloudType>::postEvolve()
{
    CloudType::postEvolve();
}


template<class CloudType>
void Foam::ElectrokineticCloud<CloudType>::autoMap(const mapPolyMesh& mapper)
{
    typedef typename particle::TrackingData<ElectrokineticCloud<CloudType> > tdType;

    tdType td(*this);

    Cloud<parcelType>::template autoMap<tdType>(td, mapper);

    this->updateMesh();
}


// ************************************************************************* //
