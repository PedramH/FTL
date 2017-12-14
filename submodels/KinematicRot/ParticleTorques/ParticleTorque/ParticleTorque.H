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

Class
    Foam::ParticleTorque

Description
    Abstract base class for particle torques

SourceFiles
    ParticleTorqueI.H
    ParticleTorque.C
    ParticleTorqueNew.C

\*---------------------------------------------------------------------------*/

#ifndef ParticleTorque_H
#define ParticleTorque_H

#include "dictionary.H"
#include "torqueSuSp.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class ParticleTorque Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class ParticleTorque
{
    // Private data

        //- Reference to the owner cloud
        CloudType& owner_;

        //- Reference to the mesh database
        const fvMesh& mesh_;

        //- Torque coefficients dictaionary
        const dictionary coeffs_;


public:

    //- Runtime type information
    TypeName("particleTorque");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        ParticleTorque,
        dictionary,
        (
            CloudType& owner,
            const fvMesh& mesh,
            const dictionary& dict
        ),
        (owner, mesh, dict)
    );


    //- Convenience typedef for return type
    typedef VectorSpace<Vector<vector>, vector, 2> returnType;


    // Constructors

        //- Construct from mesh
        ParticleTorque
        (
            CloudType& owner,
            const fvMesh& mesh,
            const dictionary& dict,
            const word& torqueType,
            const bool readCoeffs
        );

        //- Construct copy
        ParticleTorque(const ParticleTorque& pf);

        //- Construct and return a clone
        virtual autoPtr<ParticleTorque<CloudType> > clone() const
        {
            return autoPtr<ParticleTorque<CloudType> >
            (
                new ParticleTorque<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~ParticleTorque();


    //- Selector
    static autoPtr<ParticleTorque<CloudType> > New
    (
        CloudType& owner,
        const fvMesh& mesh,
        const dictionary& dict,
        const word& torqueType
    );


    // Member Functions

        // Access

            //- Return const access to the cloud owner
            inline const CloudType& owner() const;

            //- Return refernce to the cloud owner
            inline CloudType& owner();

            //- Return the mesh database
            inline const fvMesh& mesh() const;

            //- Return the torque coefficients dictionary
            inline const dictionary& coeffs() const;


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
				scalar& taopR, //added
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

            //- Return the added momentOfInertia
            virtual scalar momentOfInertiaAdd
            (
                const typename CloudType::parcelType& p,
                const scalar momentOfInertia
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ParticleTorqueI.H"

#ifdef NoRepository
    #include "ParticleTorque.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeParticleTorqueModel(CloudType)                                     \
                                                                              \
    typedef CloudType::kinematicRotCloudType kinematicRotCloudType;                 \
    defineNamedTemplateTypeNameAndDebug(ParticleTorque<kinematicRotCloudType>, 0);\
    defineTemplateRunTimeSelectionTable                                       \
    (                                                                         \
        ParticleTorque<kinematicRotCloudType>,                                    \
        dictionary                                                            \
    );


#define makeParticleTorqueModelType(SS, CloudType)                             \
                                                                              \
    typedef CloudType::kinematicRotCloudType kinematicRotCloudType;                 \
    defineNamedTemplateTypeNameAndDebug(SS<kinematicRotCloudType>, 0);           \
                                                                              \
    ParticleTorque<kinematicRotCloudType>::                                       \
        adddictionaryConstructorToTable<SS<kinematicRotCloudType> >              \
        add##SS##CloudType##kinematicRotCloudType##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
