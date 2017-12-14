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

Class
    Foam::ParticleTorqueList

Description
    List of particle torques

SourceFiles
    ParticleTorqueListI.H
    ParticleTorqueList.C

\*---------------------------------------------------------------------------*/

#ifndef ParticleTorqueList_H
#define ParticleTorqueList_H

#include "ParticleTorque.H"
#include "torqueSuSp.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class ParticleTorqueList Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class ParticleTorqueList
:
    public PtrList<ParticleTorque<CloudType> >
{
    // Private data

        //- Reference to the owner cloud
        CloudType& owner_;

        //- Reference to the mesh database
        const fvMesh& mesh_;

        //- Torques dictionary
        const dictionary dict_;

        //- Calculate coupled torques flag
        bool calcCoupled_;

        //- Calculate non-coupled torques flag
        bool calcNonCoupled_;


public:

    // Constructors

        //- Null constructor
        ParticleTorqueList(CloudType& owner, const fvMesh& mesh);

        //- Construct from mesh
        ParticleTorqueList
        (
            CloudType& owner,
            const fvMesh& mesh,
            const dictionary& dict,
            const bool readFields
        );

        //- Construct copy
        ParticleTorqueList(const ParticleTorqueList& pfl);


    //- Destructor
    virtual ~ParticleTorqueList();


    // Member Functions

        // Access

            //- Return const access to the cloud owner
            inline const CloudType& owner() const;

            //- Return refernce to the cloud owner
            inline CloudType& owner();

            //- Return the mesh database
            inline const fvMesh& mesh() const;

            //- Return the torques dictionary
            inline const dictionary& dict() const;

            //- Set the calcCoupled flag
            inline void setCalcCoupled(bool flag);

            //- Set the calcNonCoupled flag
            inline void setCalcNonCoupled(bool flag);


        // Evaluation

            //- Cache fields
            virtual void cacheFields(const bool store);

            //- Calculate the coupled torque
            virtual torqueSuSp calcCoupled
            (
                const typename CloudType::parcelType& p,
                const scalar dt,
                const scalar momentOfInertia,
                const scalar Re,
                const scalar muc,
				scalar& taop, //added
				vector& Omegac //added
            ) const;

            //- Calculate the non-coupled torque
            virtual torqueSuSp calcNonCoupled
            (
                const typename CloudType::parcelType& p,
                const scalar dt,
                const scalar momentOfInertia,
                const scalar Re,
                const scalar muc
            ) const;

            //- Return the effective momentOfInertia
            /*virtual scalar momentOfInertiaEff
            (
                const typename CloudType::parcelType& p,
                const scalar momentOfInertia
            ) const;*/
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ParticleTorqueListI.H"

#ifdef NoRepository
    #include "ParticleTorqueList.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
