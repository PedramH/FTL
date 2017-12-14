/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2013 OpenFOAM Foundation
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

#include "RubinowKellerTorque.H"
#include "mathematicalConstants.H"

using namespace Foam::constant;

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class CloudType>
Foam::scalar Foam::RubinowKellerTorque<CloudType>::RubinowKellerTorque::CR
(
    const typename CloudType::parcelType& p,
    const vector& curlUc,
    const scalar Re,
    const scalar muc
) const
{
    vector Omegad = 0.5*curlUc - p.Omega();
	scalar ReR = p.rhoc()*mag(Omegad)*sqr(p.d())/(muc + ROOTVSMALL);

	scalar CR;
	if (ReR <= 32)
	{
		CR = 64.0*mathematical::pi/(ReR+ROOTVSMALL);
	}
	else
	{
		CR = 12.9/(sqrt(ReR)+ROOTVSMALL) + 128.4/(ReR+ROOTVSMALL);
	}

    return CR;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::RubinowKellerTorque<CloudType>::RubinowKellerTorque
(
    CloudType& owner,
    const fvMesh& mesh,
    const dictionary& dict,
    const word& torqueType
)
:
    LiftTorque<CloudType>(owner, mesh, dict, torqueType)
{}


template<class CloudType>
Foam::RubinowKellerTorque<CloudType>::RubinowKellerTorque
(
	const RubinowKellerTorque& lf
)
:
    LiftTorque<CloudType>(lf)
{}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::RubinowKellerTorque<CloudType>::~RubinowKellerTorque()
{}


// ************************************************************************* //
