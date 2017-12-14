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

#include "ParticleTorqueList.H"
#include "entry.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ParticleTorqueList<CloudType>::ParticleTorqueList
(
    CloudType& owner,
    const fvMesh& mesh
)
:
    PtrList<ParticleTorque<CloudType> >(),
    owner_(owner),
    mesh_(mesh),
    dict_(dictionary::null),
    calcCoupled_(true),
    calcNonCoupled_(true)
{}


template<class CloudType>
Foam::ParticleTorqueList<CloudType>::ParticleTorqueList
(
    CloudType& owner,
    const fvMesh& mesh,
    const dictionary& dict,
    const bool readFields
)
:
    PtrList<ParticleTorque<CloudType> >(),
    owner_(owner),
    mesh_(mesh),
    dict_(dict),
    calcCoupled_(true),
    calcNonCoupled_(true)
{
    if (readFields)
    {
        wordList modelNames(dict.toc());

        Info<< "Constructing particle torques" << endl;

        if (modelNames.size() > 0)
        {
            this->setSize(modelNames.size());

            label i = 0;
            forAllConstIter(IDLList<entry>, dict, iter)
            {
                const word& model = iter().keyword();
                if (iter().isDict())
                {
                    this->set
                    (
                        i++,
                        ParticleTorque<CloudType>::New
                        (
                            owner,
                            mesh,
                            iter().dict(),
                            model
                        )
                    );
                }
                else
                {
                    this->set
                    (
                        i++,
                        ParticleTorque<CloudType>::New
                        (
                            owner,
                            mesh,
                            dict,
                            model
                        )
                    );
                }
            }
        }
        else
        {
            Info<< "    none" << endl;
        }
    }
}


template<class CloudType>
Foam::ParticleTorqueList<CloudType>::ParticleTorqueList
(
    const ParticleTorqueList& pf
)
:
    PtrList<ParticleTorque<CloudType> >(pf),
    owner_(pf.owner_),
    mesh_(pf.mesh_),
    dict_(pf.dict_)
{}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ParticleTorqueList<CloudType>::~ParticleTorqueList()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::ParticleTorqueList<CloudType>::cacheFields(const bool store)
{
    forAll(*this, i)
    {
        this->operator[](i).cacheFields(store);
    }
}


template<class CloudType>
Foam::torqueSuSp Foam::ParticleTorqueList<CloudType>::calcCoupled
(
    const typename CloudType::parcelType& p,
    const scalar dt,
    const scalar momentOfInertia,
    const scalar Re,
    const scalar muc,
	scalar& taop, //added
	vector& Omegac //added
) const
{
    torqueSuSp value(vector::zero, 0.0);

    if (calcCoupled_)
    {
        forAll(*this, i)
        {
            value += this->operator[](i).calcCoupled(p, dt, momentOfInertia, Re, muc, taop, Omegac);
        }
    }

    return value;
}


template<class CloudType>
Foam::torqueSuSp Foam::ParticleTorqueList<CloudType>::calcNonCoupled
(
    const typename CloudType::parcelType& p,
    const scalar dt,
    const scalar momentOfInertia,
    const scalar Re,
    const scalar muc
) const
{
    torqueSuSp value(vector::zero, 0.0);

    if (calcNonCoupled_)
    {
        forAll(*this, i)
        {
            value += this->operator[](i).calcNonCoupled(p, dt, momentOfInertia, Re, muc);
        }
    }

    return value;
}

/* for torques including dwp/dt
template<class CloudType>
Foam::scalar Foam::ParticleTorqueList<CloudType>::momentOfInertiaEff
(
    const typename CloudType::parcelType& p,
    const scalar momentOfInertia
) const
{
    scalar momentOfInertiaEff = momentOfInertia;
    forAll(*this, i)
    {
        momentOfInertiaEff += this->operator[](i).momentOfInertiaAdd(p, momentOfInertia);
    }

    return momentOfInertiaEff;
}
*/

// ************************************************************************* //
