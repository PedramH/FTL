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

//#include "mathematicalConstants.H"

//using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
inline Foam::ElectrokineticParcel<ParcelType>::constantProperties::constantProperties()
:
    ParcelType::constantProperties(),
    z0_(this->dict_, 0)
{}


template<class ParcelType>
inline Foam::ElectrokineticParcel<ParcelType>::constantProperties::constantProperties
(
    const constantProperties& cp
)
:
    ParcelType::constantProperties(cp),
    z0_(cp.z0_)
{}


template<class ParcelType>
inline Foam::ElectrokineticParcel<ParcelType>::constantProperties::constantProperties
(
    const dictionary& parentDict
)
:
    ParcelType::constantProperties(parentDict),
    z0_(this->dict_, "z0")
{}



template<class ParcelType>
inline Foam::ElectrokineticParcel<ParcelType>::ElectrokineticParcel
(
    const polyMesh& mesh,
    const vector& position,
    const label cellI,
    const label tetFaceI,
    const label tetPtI
)
:
    ParcelType(mesh, position, cellI, tetFaceI, tetPtI),
    z_(0),
	Ec_(vector::zero)
{}


template<class ParcelType>
inline Foam::ElectrokineticParcel<ParcelType>::ElectrokineticParcel
(
    const polyMesh& mesh,
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
    ParcelType
    (
        mesh,
        position,
        cellI,
        tetFaceI,
        tetPtI,
        typeId,
        nParticle0,
        d0,
        dTarget0,
        U0,
        constProps
    ),
    z_(constProps.z0()),
	Ec_(vector::zero)
{}


// * * * * * * * * * constantProperties Member Functions * * * * * * * * * * //

template<class ParcelType>
inline Foam::scalar
Foam::ElectrokineticParcel<ParcelType>::constantProperties::z0() const
{
    return z0_.value();
}


// * * * * * * * * * * ElectrokineticParcel Member Functions * * * * * * * * * * * * //

template<class ParcelType>
inline const Foam::scalar& Foam::ElectrokineticParcel<ParcelType>::z() const
{
    return z_;
}


template<class ParcelType>
inline Foam::scalar& Foam::ElectrokineticParcel<ParcelType>::z()
{
    return z_;
}


template<class ParcelType>
inline Foam::vector Foam::ElectrokineticParcel<ParcelType>::Ec() const
{
    return Ec_;
}


// ************************************************************************* //
