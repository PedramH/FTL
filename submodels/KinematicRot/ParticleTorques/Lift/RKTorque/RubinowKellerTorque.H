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
    Foam::RubinowKellerTorque

Description
    Magnus particle lift torque model applicable to spherical particles.

SourceFiles
    RubinowKellerTorque.C

\*---------------------------------------------------------------------------*/

#ifndef RubinowKellerTorque_H
#define RubinowKellerTorque_H

#include "LiftTorque.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class LiftTorque Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class RubinowKellerTorque
:
    public LiftTorque<CloudType>
{
protected:

    // Protected data

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

	//- Runtime type information
	TypeName("RubinowKellerTorque");

    // Constructors

        //- Construct from mesh
        RubinowKellerTorque
        (
            CloudType& owner,
            const fvMesh& mesh,
            const dictionary& dict,
            const word& torqueType = typeName
        );

        //- Construct copy
        RubinowKellerTorque(const RubinowKellerTorque& lf);

        //- Construct and return a clone
        virtual autoPtr<ParticleTorque<CloudType> > clone() const
        {
            return autoPtr<ParticleTorque<CloudType> >
            (
                new RubinowKellerTorque<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~RubinowKellerTorque();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "RubinowKellerTorque.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
