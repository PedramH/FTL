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

#include "LiftTorque.H"
#include "fvcCurl.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class CloudType>
Foam::scalar Foam::LiftTorque<CloudType>::LiftTorque::CR
(
    const typename CloudType::parcelType& p,
    const vector& curlUc,
    const scalar Re,
    const scalar muc
) const
{
    // dummy
    return 0.0;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::LiftTorque<CloudType>::LiftTorque
(
    CloudType& owner,
    const fvMesh& mesh,
    const dictionary& dict,
    const word& torqueType
)
:
    ParticleTorque<CloudType>(owner, mesh, dict, torqueType, true),
    UName_(this->coeffs().template lookupOrDefault<word>("U", "U")),
    curlUcInterpPtr_(NULL)
{}


template<class CloudType>
Foam::LiftTorque<CloudType>::LiftTorque(const LiftTorque& lf)
:
    ParticleTorque<CloudType>(lf),
    UName_(lf.UName_),
    curlUcInterpPtr_(NULL)
{}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::LiftTorque<CloudType>::~LiftTorque()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::LiftTorque<CloudType>::cacheFields(const bool store)
{
    static word fName("curlUcDt");

    bool fieldExists = this->mesh().template foundObject<volVectorField>(fName);

    if (store)
    {
        if (!fieldExists)
        {
            const volVectorField& Uc = this->mesh().template
                lookupObject<volVectorField>(UName_);

            volVectorField* curlUcPtr =
                new volVectorField(fName, fvc::curl(Uc));

            curlUcPtr->store();
        }

        const volVectorField& curlUc = this->mesh().template
            lookupObject<volVectorField>(fName);

        curlUcInterpPtr_.reset
        (
            interpolation<vector>::New
            (
                this->owner().solution().interpolationSchemes(),
                curlUc
            ).ptr()
        );
    }
    else
    {
        curlUcInterpPtr_.clear();

        if (fieldExists)
        {
            const volVectorField& curlUc = this->mesh().template
                lookupObject<volVectorField>(fName);

            const_cast<volVectorField&>(curlUc).checkOut();
        }
    }
}


template<class CloudType>
Foam::torqueSuSp Foam::LiftTorque<CloudType>::calcCoupled
(
    const typename CloudType::parcelType& p,
    const scalar dt,
    const scalar momentOfInertia,
    const scalar Re,
    const scalar muc,
	scalar& taopR, //added
	vector& Omegac //added
) const
{
    torqueSuSp value(vector::zero, 0.0);

    vector curlUc =
        curlUcInterp().interpolate(p.position(), p.currentTetIndices());

    Omegac = 0.5*curlUc;

    scalar CR = this->CR(p, curlUc, Re, muc);

	vector Omegad = Omegac - p.Omega();
	scalar Sp = p.rho()/2*pow(p.d()/2,5)*CR*mag(Omegad);

	taopR = momentOfInertia / (Sp+ROOTVSMALL);

    value.Sp() = Sp;

    return value;
}


// ************************************************************************* //
