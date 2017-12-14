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
    Foam::MagnusLiftRotForce

Description
    Magnus particle lift force model applicable to spherical particles.

SourceFiles
    MagnusLiftRotForce.C

\*---------------------------------------------------------------------------*/

#ifndef MagnusLiftRotForce_H
#define MagnusLiftRotForce_H

#include "LiftRotForce.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class MagnusLiftRotForce Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class MagnusLiftRotForce
:
    public LiftRotForce<CloudType>
{
protected:

    // Protected Member Functions

        //- Calculate the lift coefficient
        virtual scalar ClRA
        (
            const typename CloudType::parcelType& p,
            const vector& curlUc,
            const scalar Re,
            const scalar muc
        ) const;


public:

    //- Runtime type information
    TypeName("MagnusLiftRotForce");


    // Constructors

        //- Construct from mesh
        MagnusLiftRotForce
        (
            CloudType& owner,
            const fvMesh& mesh,
            const dictionary& dict,
            const word& forceType = typeName
        );

        //- Construct copy
        MagnusLiftRotForce(const MagnusLiftRotForce& lf);

        //- Construct and return a clone
        virtual autoPtr<ParticleForce<CloudType> > clone() const
        {
            return autoPtr<ParticleForce<CloudType> >
            (
                new MagnusLiftRotForce<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~MagnusLiftRotForce();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "MagnusLiftRotForce.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
