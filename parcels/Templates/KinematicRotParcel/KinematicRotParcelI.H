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
inline Foam::KinematicRotParcel<ParcelType>::constantProperties::constantProperties()
:
    ParcelType::constantProperties(),
    Omega0_(this->dict_, vector::zero)
{}


template<class ParcelType>
inline Foam::KinematicRotParcel<ParcelType>::constantProperties::constantProperties
(
    const constantProperties& cp
)
:
    ParcelType::constantProperties(cp),
    Omega0_(cp.Omega0_)
{}


template<class ParcelType>
inline Foam::KinematicRotParcel<ParcelType>::constantProperties::constantProperties
(
    const dictionary& parentDict
)
:
    ParcelType::constantProperties(parentDict),
    Omega0_(this->dict_, "Omega0")
{}


template<class ParcelType>
inline Foam::KinematicRotParcel<ParcelType>::KinematicRotParcel
(
    const polyMesh& mesh,
    const vector& position,
    const label cellI,
    const label tetFaceI,
    const label tetPtI
)
:
    ParcelType(mesh, position, cellI, tetFaceI, tetPtI),
    Omega_(vector::zero),
	Omegac_(vector::zero)
{}


template<class ParcelType>
inline Foam::KinematicRotParcel<ParcelType>::KinematicRotParcel
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
    Omega_(constProps.Omega0()),
	Omegac_(vector::zero)
{}


// * * * * * * * * * constantProperties Member Functions * * * * * * * * * * //

template<class ParcelType>
inline Foam::vector
Foam::KinematicRotParcel<ParcelType>::constantProperties::Omega0() const
{
    return Omega0_.value();
}


// * * * * * * * * * * KinematicRotParcel Member Functions * * * * * * * * * * * * //

template<class ParcelType>
inline const Foam::vector& Foam::KinematicRotParcel<ParcelType>::Omega() const
{
    return Omega_;
}


template<class ParcelType>
inline Foam::vector& Foam::KinematicRotParcel<ParcelType>::Omega()
{
    return Omega_;
}


template<class ParcelType>
inline Foam::vector Foam::KinematicRotParcel<ParcelType>::Omegac() const
{
    return Omegac_;
}


template<class ParcelType>
inline Foam::scalar Foam::KinematicRotParcel<ParcelType>::momentOfInertiaCell 
(
    const label cellI
) const
{
	scalar coef = 1.0; //correcion coefficient, increasing his coefficient reduces the Omegac change during integration

	//sphere (minimum probable)
	//scalar dE = pow(6/pi*this->mesh().cellVolumes()[cellI],1.0/3);
    //scalar MOI = 0.1*this->massCell(cellI)*sqr(dE);
	//cube
	scalar aE = pow(this->mesh().cellVolumes()[cellI],1.0/3);
    scalar MOI = 1.0/6*this->massCell(cellI)*sqr(aE);

	return coef*MOI;
}


// ************************************************************************* //
