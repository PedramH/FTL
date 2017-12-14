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
    Foam::kinematicRotCloud

Description
    Virtual abstract base class for templated ThermoCloud

SourceFiles
    kinematicRotCloud.C

\*---------------------------------------------------------------------------*/

#ifndef kinematicRotCloud_H
#define kinematicRotCloud_H

#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class kinematicRotCloud Declaration
\*---------------------------------------------------------------------------*/

class kinematicRotCloud
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        kinematicRotCloud(const kinematicRotCloud&);

        //- Disallow default bitwise assignment
        void operator=(const kinematicRotCloud&);


public:

    //- Runtime type information
    TypeName("kinematicRotCloud");

    // Constructors

        //- Null constructor
        kinematicRotCloud();


    //- Destructor
    virtual ~kinematicRotCloud();


    // Member Functions

        // Check

            // Radiation

                //- Total rotational kinetic energy in the system
//            virtual scalar rotationalKineticEnergyOfSystem() const = 0;

            //- Penetration for percentage of the current total mass
//            virtual scalar penetration(const scalar& fraction) const = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
