/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::torqueSuSp

Description
    Helper container for torque Su and Sp terms.

        F = Sp(U - Up) + Su

    Explicit contribution, Su specified as a torque
    Implicit coefficient, Sp specified as torque/velocity

SourceFiles
    torqueSuSpI.H

\*---------------------------------------------------------------------------*/

#ifndef torqueSuSp_H
#define torqueSuSp_H

#include "Tuple2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions
class torqueSuSp;
inline torqueSuSp operator+(const torqueSuSp& susp1, const torqueSuSp& susp2);
inline torqueSuSp operator*(const torqueSuSp& susp, const scalar s);

/*---------------------------------------------------------------------------*\
                          Class torqueSuSp Declaration
\*---------------------------------------------------------------------------*/

class torqueSuSp
:
    public Tuple2<vector, scalar>
{

public:

    // Constructors

        //- Construct null
        inline torqueSuSp();

        //- Construct given Tuple2
        inline torqueSuSp(const Tuple2<vector, scalar>& susp);

        //- Construct given two components
        inline torqueSuSp(const vector& Su, const scalar Sp);

        //- Construct from Istream
        inline torqueSuSp(Istream&);


    // Member Functions

        // Access

            //- Return const access to the explicit contribution [kg.m/s2]
            inline const vector& Su() const;

            //- Return const access to the implicit coefficient [kg/s]
            inline scalar Sp() const;


        // Edit

            //- Return reference to the explicit contribution
            inline vector& Su();

            //- Return reference to the implicit coefficient
            inline scalar& Sp();


        // Operators

            //- Assignment
            inline void operator=(const torqueSuSp& susp);

            //- Addition
            inline void operator+=(const torqueSuSp& susp);

            //- Subtraction
            inline void operator-=(const torqueSuSp& susp);


        // Friend operators

            //- Addition
            friend inline torqueSuSp operator*
            (
                const torqueSuSp& susp1,
                const torqueSuSp& susp2
            );

            //- Multiplication
            friend inline torqueSuSp operator*
            (
                const scalar s,
                const torqueSuSp& susp
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "torqueSuSpI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
