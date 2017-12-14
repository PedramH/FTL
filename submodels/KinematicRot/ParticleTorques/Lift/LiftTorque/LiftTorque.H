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

Class
    Foam::LiftTorque

Description
    Base class for particle lift torque models

SourceFiles
    LiftTorqueI.H
    LiftTorque.C

\*---------------------------------------------------------------------------*/

#ifndef LiftTorque_H
#define LiftTorque_H

#include "ParticleTorque.H"
#include "volFields.H"
#include "interpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class LiftTorque Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class LiftTorque
:
    public ParticleTorque<CloudType>
{
protected:

    // Protected data

        //- Name of velocity field
        const word UName_;

        //- Curk of carrier phase velocity interpolator
        autoPtr<interpolation<vector> > curlUcInterpPtr_;


    // Potected Member Functions

        //- Calculate the rotational lift coefficient
        virtual scalar CR
        (
            const typename CloudType::parcelType& p,
            const vector& curlUc,
            const scalar Re,
            const scalar muc
        ) const;


public:

    // Constructors

        //- Construct from mesh
        LiftTorque
        (
            CloudType& owner,
            const fvMesh& mesh,
            const dictionary& dict,
            const word& torqueType
        );

        //- Construct copy
        LiftTorque(const LiftTorque& lf);

        //- Construct and return a clone
        virtual autoPtr<ParticleTorque<CloudType> > clone() const
        {
            return autoPtr<ParticleTorque<CloudType> >
            (
                new LiftTorque<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~LiftTorque();


    // Member Functions

        // Access

            //- Return the curl of the carrier phase velocity interpolator
            inline const interpolation<vector>& curlUcInterp() const;


        // Evaluation

            //- Cache fields
            virtual void cacheFields(const bool store);

            //- Calculate the non-coupled torque
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
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "LiftTorqueI.H"

#ifdef NoRepository
    #include "LiftTorque.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
