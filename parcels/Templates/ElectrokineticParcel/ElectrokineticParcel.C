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

#include "ElectrokineticParcel.H"

// * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * * //

template<class ParcelType>
template<class TrackData>
void Foam::ElectrokineticParcel<ParcelType>::setCellValues
(
    TrackData& td,
    const scalar dt,
    const label cellI
)
{
    ParcelType::setCellValues(td, dt, cellI);

    tetIndices tetIs = this->currentTetIndices(); 

    Ec_ = td.EInterp().interpolate(this->position(), tetIs); 
}


template<class ParcelType>
template<class TrackData>
void Foam::ElectrokineticParcel<ParcelType>::cellValueSourceCorrection
(
    TrackData& td,
    const scalar dt,
    const label cellI
)
{
	ParcelType::cellValueSourceCorrection(td, dt, cellI);
	//if the change in Ec_ can be approximated from ETrans, the calculation will be improved.
    //Ec_ += td.cloud().ETrans()[cellI]/????(cellI); 
}


template<class ParcelType>
template<class TrackData>
void Foam::ElectrokineticParcel<ParcelType>::calc
(
    TrackData& td,
    const scalar dt,
    const label cellI
)
{
	ParcelType::calc(td, dt, cellI);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::ElectrokineticParcel<ParcelType>::ElectrokineticParcel
(
    const ElectrokineticParcel<ParcelType>& p
)
:
    ParcelType(p),
    z_(p.z_),
	Ec_(p.Ec_)
{}


template<class ParcelType>
Foam::ElectrokineticParcel<ParcelType>::ElectrokineticParcel
(
    const ElectrokineticParcel<ParcelType>& p,
    const polyMesh& mesh
)
:
    ParcelType(p, mesh),
    z_(p.z_),
	Ec_(p.Ec_)
{}


// * * * * * * * * * * * * * * IOStream operators  * * * * * * * * * * * * * //

#include "ElectrokineticParcelIO.C"

// ************************************************************************* //
