/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "HybridModulation.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::HybridModulation<CloudType>::HybridModulation
(
    const dictionary& dict,
    CloudType& owner
)
:
    ModulationModel<CloudType>(dict, owner, typeName), //ModulationModel<CloudType>(owner),
	alpha_(this->coeffDict().lookupOrDefault("alpha",1.0)),
	beta_(this->coeffDict().lookupOrDefault("beta",0.5)),
	UpImplicit_(this->coeffDict().lookupOrDefault("UpImplicit",false))
{
	//change solution flag to compute and use turbulence modulation
	if (owner.solution().coupled()) owner.solution().turbulenceCoupling() = true;
}


template<class CloudType>
Foam::HybridModulation<CloudType>::HybridModulation
(
    HybridModulation<CloudType>& dm
)
:
    ModulationModel<CloudType>(dm),
	alpha_(dm.alpha_),
	beta_(dm.beta_)
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::HybridModulation<CloudType>::~HybridModulation()
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::symmTensor Foam::HybridModulation<CloudType>::update
(
    const vector UcMean,
	const vector Uc,
	const vector Up0,
	const vector Upn,
	const vector Su
)
{
	vector Up = Up0; 
	if (UpImplicit_) Up = Upn;

	vector uIFluc = beta_*(Up - UcMean) + alpha_*(Uc - UcMean);

	// Reynolds stress transfer from the particle to the carrier phase
	return 2 * symm(uIFluc * Su);
}


// ************************************************************************* //
