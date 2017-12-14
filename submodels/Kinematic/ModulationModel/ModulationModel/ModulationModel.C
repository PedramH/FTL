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

#include "ModulationModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ModulationModel<CloudType>::ModulationModel(CloudType& owner)
:
    CloudSubModelBase<CloudType>(owner)
{}


template<class CloudType>
Foam::ModulationModel<CloudType>::ModulationModel
(
    const dictionary& dict,
    CloudType& owner,
    const word& type
)
:
    CloudSubModelBase<CloudType>(owner, dict, typeName, type)
{}


template<class CloudType>
Foam::ModulationModel<CloudType>::ModulationModel
(
    ModulationModel<CloudType>& dm
)
:
    CloudSubModelBase<CloudType>(dm)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ModulationModel<CloudType>::~ModulationModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::symmTensor Foam::ModulationModel<CloudType>::update
(
    const vector UcMean,
	const vector Uc,
	const vector Up0,
	const vector Upn,
	const vector Su
)
{
    notImplemented
    (
        "Foam::symmTensor Foam::ModulationModel<CloudType>::update"
        "("
            "const vector, "
            "const vector,"
			"const vector,"
			"const vector,"
			"const vector"
        ")"
    );

    return symmTensor::zero;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ModulationModelNew.C"

// ************************************************************************* //
