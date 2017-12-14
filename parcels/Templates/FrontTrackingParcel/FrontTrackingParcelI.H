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

#include "mathematicalConstants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
inline
Foam::FrontTrackingParcel<ParcelType>::constantProperties::constantProperties()
:
    dict_(dictionary::null),
    parcelTypeId_(dict_, -1),
    rhoMin_(dict_, 0.0),
    rho0_(dict_, 0.0),
    minParticleMass_(dict_, 0.0)
{}


template<class ParcelType>
inline Foam::FrontTrackingParcel<ParcelType>::constantProperties::constantProperties
(
    const constantProperties& cp
)
:
    dict_(cp.dict_),
    parcelTypeId_(cp.parcelTypeId_),
    rhoMin_(cp.rhoMin_),
    rho0_(cp.rho0_),
    minParticleMass_(cp.minParticleMass_)
{}


template<class ParcelType>
inline Foam::FrontTrackingParcel<ParcelType>::constantProperties::constantProperties
(
    const dictionary& parentDict
)
:
    dict_(parentDict.subOrEmptyDict("constantProperties")),
    parcelTypeId_(dict_, "parcelTypeId", -1),
    rhoMin_(dict_, "rhoMin", 1e-15),
    rho0_(dict_, "rho0"),
    minParticleMass_(dict_, "minParticleMass", 1e-15)
{}


template<class ParcelType>
inline Foam::FrontTrackingParcel<ParcelType>::FrontTrackingParcel
(
    const polyMesh& owner,
    const vector& position,
    const label cellI,
    const label tetFaceI,
    const label tetPtI
)
:
    ParcelType(owner, position, cellI, tetFaceI, tetPtI),
    active_(true),
    typeId_(-1),
    nParticle_(0),
    d_(0.0),
    dTarget_(0.0),
    U_(vector::zero),
    rho_(0.0),
    age_(0.0),
    tTurb_(0.0),
    UTurb_(vector::zero),
    currentIndex_(0),//Jafari added
    bubbleIndex_(0),//Jafari added
    rhoc_(0.0),
    Uc_(vector::zero),
	UcMean_(vector::zero), //added ttw
    muc_(0.0)
{}


template<class ParcelType>
inline Foam::FrontTrackingParcel<ParcelType>::FrontTrackingParcel
(
    const polyMesh& owner,
    const vector& position,
    const label cellI,
    const label tetFaceI,
    const label tetPtI,
    const label typeId,
    const scalar nParticle0,
    const scalar d0,
    const scalar dTarget0,
    const vector& U0,
    const constantProperties& constProps
)
:
    ParcelType(owner, position, cellI, tetFaceI, tetPtI),
    active_(true),
    typeId_(typeId),
    nParticle_(nParticle0),
    d_(d0),
    dTarget_(dTarget0),
    U_(U0),
    rho_(constProps.rho0()),
    age_(0.0),
    tTurb_(0.0),
    UTurb_(vector::zero),
    currentIndex_(0),//Jafari added
    bubbleIndex_(0),//Jafari added
    rhoc_(0.0),
    Uc_(vector::zero),
	UcMean_(vector::zero), //added ttw
    muc_(0.0)
{}


// * * * * * * * * * constantProperties Member Functions * * * * * * * * * * //

template<class ParcelType>
inline const Foam::dictionary&
Foam::FrontTrackingParcel<ParcelType>::constantProperties::dict() const
{
    return dict_;
}


template<class ParcelType>
inline Foam::label
Foam::FrontTrackingParcel<ParcelType>::constantProperties::parcelTypeId() const
{
    return parcelTypeId_.value();
}


template<class ParcelType>
inline Foam::scalar
Foam::FrontTrackingParcel<ParcelType>::constantProperties::rhoMin() const
{
    return rhoMin_.value();
}


template<class ParcelType>
inline Foam::scalar
Foam::FrontTrackingParcel<ParcelType>::constantProperties::rho0() const
{
    return rho0_.value();
}


template<class ParcelType>
inline Foam::scalar
Foam::FrontTrackingParcel<ParcelType>::constantProperties::minParticleMass() const
{
    return minParticleMass_.value();
}


// * * * * * * * FrontTrackingParcel Member Functions  * * * * * * * //

template<class ParcelType>
inline bool Foam::FrontTrackingParcel<ParcelType>::active() const
{
    return active_;
}


template<class ParcelType>
inline Foam::label Foam::FrontTrackingParcel<ParcelType>::typeId() const
{
    return typeId_;
}


template<class ParcelType>
inline Foam::scalar Foam::FrontTrackingParcel<ParcelType>::nParticle() const
{
    return nParticle_;
}


template<class ParcelType>
inline Foam::scalar Foam::FrontTrackingParcel<ParcelType>::d() const
{
    return d_;
}


template<class ParcelType>
inline Foam::scalar Foam::FrontTrackingParcel<ParcelType>::dTarget() const
{
    return dTarget_;
}


template<class ParcelType>
inline const Foam::vector& Foam::FrontTrackingParcel<ParcelType>::U() const
{
    return U_;
}


template<class ParcelType>
inline Foam::scalar Foam::FrontTrackingParcel<ParcelType>::rho() const
{
    return rho_;
}


template<class ParcelType>
inline Foam::scalar Foam::FrontTrackingParcel<ParcelType>::age() const
{
    return age_;
}


template<class ParcelType>
inline Foam::scalar Foam::FrontTrackingParcel<ParcelType>::tTurb() const
{
    return tTurb_;
}


template<class ParcelType>
inline const Foam::vector& Foam::FrontTrackingParcel<ParcelType>::UTurb() const
{
    return UTurb_;
}

//Jafari added
template<class ParcelType>
inline const Foam::label& Foam::FrontTrackingParcel<ParcelType>::currentIndex() const
{
    return currentIndex_;
}

//Jafari added
template<class ParcelType>
inline const Foam::label& Foam::FrontTrackingParcel<ParcelType>::bubbleIndex() const
{
    return bubbleIndex_;
}


template<class ParcelType>
inline Foam::scalar Foam::FrontTrackingParcel<ParcelType>::rhoc() const
{
    return rhoc_;
}


template<class ParcelType>
inline const Foam::vector& Foam::FrontTrackingParcel<ParcelType>::Uc() const
{
    return Uc_;
}


template<class ParcelType> //added ttw
inline const Foam::vector& Foam::FrontTrackingParcel<ParcelType>::UcMean() const
{
    return UcMean_;
}


template<class ParcelType>
inline Foam::scalar Foam::FrontTrackingParcel<ParcelType>::muc() const
{
    return muc_;
}


template<class ParcelType>
inline bool& Foam::FrontTrackingParcel<ParcelType>::active()
{
    return active_;
}


template<class ParcelType>
inline Foam::label& Foam::FrontTrackingParcel<ParcelType>::typeId()
{
    return typeId_;
}


template<class ParcelType>
inline Foam::scalar& Foam::FrontTrackingParcel<ParcelType>::nParticle()
{
    return nParticle_;
}


template<class ParcelType>
inline Foam::scalar& Foam::FrontTrackingParcel<ParcelType>::d()
{
    return d_;
}


template<class ParcelType>
inline Foam::scalar& Foam::FrontTrackingParcel<ParcelType>::dTarget()
{
    return dTarget_;
}


template<class ParcelType>
inline Foam::vector& Foam::FrontTrackingParcel<ParcelType>::U()
{
    return U_;
}


template<class ParcelType>
inline Foam::scalar& Foam::FrontTrackingParcel<ParcelType>::rho()
{
    return rho_;
}


template<class ParcelType>
inline Foam::scalar& Foam::FrontTrackingParcel<ParcelType>::age()
{
    return age_;
}


template<class ParcelType>
inline Foam::scalar& Foam::FrontTrackingParcel<ParcelType>::tTurb()
{
    return tTurb_;
}


template<class ParcelType>
inline Foam::vector& Foam::FrontTrackingParcel<ParcelType>::UTurb()
{
    return UTurb_;
}

//Jafari added
template<class ParcelType>
inline Foam::label& Foam::FrontTrackingParcel<ParcelType>::currentIndex()
{
    return currentIndex_;
}

//Jafari added
template<class ParcelType>
inline Foam::label& Foam::FrontTrackingParcel<ParcelType>::bubbleIndex()
{
    return bubbleIndex_;
}

template<class ParcelType>
inline Foam::label Foam::FrontTrackingParcel<ParcelType>::faceInterpolation() const
{
    // Use volume-based interpolation if dealing with external faces
    if (this->cloud().internalFace(this->face()))
    {
        return this->face();
    }
    else
    {
        return -1;
    }
}


template<class ParcelType>
inline Foam::scalar Foam::FrontTrackingParcel<ParcelType>::massCell
(
    const label cellI
) const
{
    return rhoc_*this->mesh().cellVolumes()[cellI];
}


template<class ParcelType>
inline Foam::scalar Foam::FrontTrackingParcel<ParcelType>::mass() const
{
    return rho_*volume();
}


template<class ParcelType>
inline Foam::scalar Foam::FrontTrackingParcel<ParcelType>::momentOfInertia() const
{
    return 0.1*mass()*sqr(d_);
}


template<class ParcelType>
inline Foam::scalar Foam::FrontTrackingParcel<ParcelType>::volume() const
{
    return volume(d_);
}


template<class ParcelType>
inline Foam::scalar Foam::FrontTrackingParcel<ParcelType>::volume(const scalar d)
{
    return pi/6.0*pow3(d);
}


template<class ParcelType>
inline Foam::scalar Foam::FrontTrackingParcel<ParcelType>::areaP() const
{
    return areaP(d_);
}


template<class ParcelType>
inline Foam::scalar Foam::FrontTrackingParcel<ParcelType>::areaP(const scalar d)
{
    return 0.25*areaS(d);
}


template<class ParcelType>
inline Foam::scalar Foam::FrontTrackingParcel<ParcelType>::areaS() const
{
    return areaS(d_);
}


template<class ParcelType>
inline Foam::scalar Foam::FrontTrackingParcel<ParcelType>::areaS(const scalar d)
{
    return pi*d*d;
}


template<class ParcelType>
inline Foam::scalar Foam::FrontTrackingParcel<ParcelType>::Re
(
    const vector& U,
    const scalar d,
    const scalar rhoc,
    const scalar muc
) const
{
    return rhoc*mag(U - Uc_)*d/(muc + ROOTVSMALL);
}


template<class ParcelType>
inline Foam::scalar Foam::FrontTrackingParcel<ParcelType>::We
(
    const vector& U,
    const scalar d,
    const scalar rhoc,
    const scalar sigma
) const
{
    return rhoc*magSqr(U - Uc_)*d/(sigma + ROOTVSMALL);
}


template<class ParcelType>
inline Foam::scalar Foam::FrontTrackingParcel<ParcelType>::Eo
(
    const vector& a,
    const scalar d,
    const scalar sigma
) const
{
    vector dir = U_/(mag(U_) + ROOTVSMALL);
    return mag(a & dir)*(rho_ - rhoc_)*sqr(d)/(sigma + ROOTVSMALL);
}


// ************************************************************************* //
