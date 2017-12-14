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

#include "MyKinematicLookupTableInjection.H"
#include "scalarIOList.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::MyKinematicLookupTableInjection<CloudType>::MyKinematicLookupTableInjection
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName
)
:
    InjectionModel<CloudType>(dict, owner, modelName, typeName),
	ranGen_(label(0)),
    inputFileName_(this->coeffDict().lookup("inputFile")),
    duration_(readScalar(this->coeffDict().lookup("duration"))),
    parcelsPerSecond_
    (
        readScalar(this->coeffDict().lookup("parcelsPerSecond"))
    ),
    randomise_(readBool(this->coeffDict().lookup("randomise"))),
    injectors_
    (
        IOobject
        (
            inputFileName_,
            owner.db().time().constant(),
            owner.db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    injectorCells_(0),
    injectorTetFaces_(0),
    injectorTetPts_(0)
{
    duration_ = owner.db().time().userTimeToTime(duration_);

    // Set/cache the injector cells
    injectorCells_.setSize(injectors_.size());
    injectorTetFaces_.setSize(injectors_.size());
    injectorTetPts_.setSize(injectors_.size());

    updateMesh();

    // Determine volume of particles to inject
    this->volumeTotal_ = 0.0;
    forAll(injectors_, i)
    {
        this->volumeTotal_ += injectors_[i].mDot()/injectors_[i].rho();
    }
    this->volumeTotal_ *= duration_;
}


template<class CloudType>
Foam::MyKinematicLookupTableInjection<CloudType>::MyKinematicLookupTableInjection
(
    const MyKinematicLookupTableInjection<CloudType>& im
)
:
    InjectionModel<CloudType>(im),
	ranGen_(im.ranGen_),
    inputFileName_(im.inputFileName_),
    duration_(im.duration_),
    parcelsPerSecond_(im.parcelsPerSecond_),
    randomise_(im.randomise_),
    injectors_(im.injectors_),
    injectorCells_(im.injectorCells_),
    injectorTetFaces_(im.injectorTetFaces_),
    injectorTetPts_(im.injectorTetPts_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::MyKinematicLookupTableInjection<CloudType>::~MyKinematicLookupTableInjection()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::MyKinematicLookupTableInjection<CloudType>::updateMesh()
{
    // Set/cache the injector cells
    forAll(injectors_, i)
    {
        this->findCellAtPosition
        (
            injectorCells_[i],
            injectorTetFaces_[i],
            injectorTetPts_[i],
            injectors_[i].x()
        );
    }
}


template<class CloudType>
Foam::scalar Foam::MyKinematicLookupTableInjection<CloudType>::timeEnd() const
{
    return this->SOI_ + duration_;
}


template<class CloudType>
Foam::label Foam::MyKinematicLookupTableInjection<CloudType>::parcelsToInject
(
    const scalar time0,
    const scalar time1
)
{
    if ((time0 >= 0.0) && (time0 < duration_))
    {
        return floor(injectorCells_.size()*(time1 - time0)*parcelsPerSecond_);
    }
    else
    {
        return 0;
    }
}


template<class CloudType>
Foam::scalar Foam::MyKinematicLookupTableInjection<CloudType>::volumeToInject
(
    const scalar time0,
    const scalar time1
)
{
    scalar volume = 0.0;
    if ((time0 >= 0.0) && (time0 < duration_))
    {
        forAll(injectors_, i)
        {
            volume += injectors_[i].mDot()/injectors_[i].rho()*(time1 - time0);
        }
    }

    return volume;
}


template<class CloudType>
void Foam::MyKinematicLookupTableInjection<CloudType>::setPositionAndCell
(
    const label parcelI,
    const label nParcels,
    const scalar time,
    vector& position,
    label& cellOwner,
    label& tetFaceI,
    label& tetPtI
)
{
    label injectorI = 0;
    if (randomise_)
    {
        cachedRandom& rnd = this->owner().rndGen();
        injectorI = rnd.position<label>(0, injectorCells_.size() - 1);
    }
    else
    {
        injectorI = parcelI*injectorCells_.size()/nParcels;
    }

    position = injectors_[injectorI].x();
    cellOwner = injectorCells_[injectorI];
    tetFaceI = injectorTetFaces_[injectorI];
    tetPtI = injectorTetPts_[injectorI];
}


template<class CloudType>
void Foam::MyKinematicLookupTableInjection<CloudType>::setProperties
(
    const label parcelI,
    const label nParcels,
    const scalar,
    typename CloudType::parcelType& parcel
)
{
    label injectorI = parcelI*injectorCells_.size()/nParcels;

	//cachedRandom& rnd = this->owner().rndGen();

    // set particle diameter
	scalar dmin = injectors_[injectorI].d() - injectors_[injectorI].deltaD()/2;

	//scalar dmax = dmin + injectors_[injectorI].deltaD();
    //parcel.d() = rnd.position<scalar>(dmin, dmax);

	parcel.d() = dmin + ranGen_.scalar01() * injectors_[injectorI].deltaD();

    // set particle density
    parcel.rho() = injectors_[injectorI].rho();

	//calculate particle velocity
	
	vector randv;
	randv.x() = ranGen_.GaussNormal();
	randv.y() = ranGen_.GaussNormal();
	randv.z() = ranGen_.GaussNormal();

	//symmTensor sisj = cmptMultiply(1 - ai * ai, injectors_[injectorI].uiuj());
	symmTensor sisj = injectors_[injectorI].uiuj();
	tensor bij = CholskyDecompose(sisj);

	vector up;
	//up = cmptMultiply(ai,up) + bij & randv;
	up.x() = bij.xx()*randv.x() + bij.xy()*randv.y() + bij.xz()*randv.z();
	up.y() = bij.yx()*randv.x() + bij.yy()*randv.y() + bij.yz()*randv.z();
	up.z() = bij.zx()*randv.x() + bij.zy()*randv.y() + bij.zz()*randv.z();

	// set particle velocity
    parcel.U() = injectors_[injectorI].U() + up;

	// calculate fluid phase flactuating velocity (isotropic)
	randv.x() = ranGen_.GaussNormal();
	randv.y() = ranGen_.GaussNormal();
	randv.z() = ranGen_.GaussNormal();

	vector uf = sqrt(2*injectors_[injectorI].kf()/3) * randv;

	// set fluid phase flactuating velocity
    parcel.UTurb() = uf;
}

template<class CloudType>
Foam::tensor Foam::MyKinematicLookupTableInjection<CloudType>::CholskyDecompose(symmTensor& sisj) 
{
	tensor bij = tensor::zero;

	bij.zz() = sqrt(sisj.zz());
	if (bij.zz() > 0) 
	{
		bij.yz() = sisj.yz() / bij.zz();
		bij.xz() = sisj.xz() / bij.zz();
	}

	bij.yy() = sqrt(mag( sisj.yy() - bij.yz()*bij.yz() ));
	if (bij.yy() > 0) bij.xy() = ( sisj.xy() - bij.yz()*bij.xz() ) / bij.yy();

	bij.xx() = sqrt(mag( sisj.xx() - bij.xz()*bij.xz() - bij.xy()*bij.xy() ));

	return bij;
}


template<class CloudType>
bool Foam::MyKinematicLookupTableInjection<CloudType>::fullyDescribed() const
{
    return true;
}


template<class CloudType>
bool Foam::MyKinematicLookupTableInjection<CloudType>::validInjection(const label)
{
    return true;
}


// ************************************************************************* //
