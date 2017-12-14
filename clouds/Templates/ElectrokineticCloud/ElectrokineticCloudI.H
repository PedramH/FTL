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

//#include "electromagneticConstants.H"

//using namespace Foam::constant;

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
inline const Foam::ElectrokineticCloud<CloudType>&
Foam::ElectrokineticCloud<CloudType>::cloudCopy() const
{
    return cloudCopyPtr_();
}


template<class CloudType>
inline const typename CloudType::particleType::constantProperties&
Foam::ElectrokineticCloud<CloudType>::constProps() const
{
    return constProps_;
}


template<class CloudType>
inline typename CloudType::particleType::constantProperties&
Foam::ElectrokineticCloud<CloudType>::constProps()
{
    return constProps_;
}


template<class CloudType>
inline const Foam::volVectorField& Foam::ElectrokineticCloud<CloudType>::E() const
{
    return E_;
}


template<class CloudType>
inline const Foam::tmp<Foam::volScalarField>
Foam::ElectrokineticCloud<CloudType>::SPhi() const
{
    tmp<volScalarField> tSPhi
    (
        new volScalarField
        (
            IOobject
            (
                this->name() + ":SPhi",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            dimensionedScalar("zero", dimensionSet(0, -3, 1, 0, 0, 1, 0), 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    volScalarField& SPhi = tSPhi();
    forAllConstIter(typename ElectrokineticCloud<CloudType>, *this, iter)
    {
        const parcelType& p = iter();
        const label cellI = p.cell();

        SPhi[cellI] += p.nParticle()*p.z();
    }

    SPhi.internalField() *= 1.60217657e-19/this->mesh().V();
    SPhi.correctBoundaryConditions();

    return tSPhi;
}


// ************************************************************************* //
