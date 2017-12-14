/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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

#include "MagnusLiftRotForce.H"
#include "mathematicalConstants.H"

using namespace Foam::constant;

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class CloudType>
Foam::scalar Foam::MagnusLiftRotForce<CloudType>::MagnusLiftRotForce::ClRA
(
    const typename CloudType::parcelType& p,
    const vector& curlUc,
    const scalar Re,
    const scalar muc
) const
{
	vector Omegad = 0.5*curlUc - p.Omega();
    scalar Rer = p.rhoc()*mag(Omegad)*sqr(p.d())/(muc + ROOTVSMALL);

    scalar ClR = 0.45 + (Rer/(Re + ROOTVSMALL) - 0.45)*exp(-0.05684*pow(Rer,0.4)*pow(Re,0.3));

    return (mathematical::pi/4*sqr(p.d()))*ClR;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::MagnusLiftRotForce<CloudType>::MagnusLiftRotForce
(
    CloudType& owner,
    const fvMesh& mesh,
    const dictionary& dict,
    const word& forceType
)
:
    LiftRotForce<CloudType>(owner, mesh, dict, forceType)
{}


template<class CloudType>
Foam::MagnusLiftRotForce<CloudType>::MagnusLiftRotForce
(
    const MagnusLiftRotForce& lf
)
:
    LiftRotForce<CloudType>(lf)
{}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::MagnusLiftRotForce<CloudType>::~MagnusLiftRotForce()
{}


// ************************************************************************* //
