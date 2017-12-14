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
    Foam::electrokineticCloud

Description
    Virtual abstract base class for templated ThermoCloud

SourceFiles
    electrokineticCloud.C

\*---------------------------------------------------------------------------*/

#ifndef electrokineticCloud_H
#define electrokineticCloud_H

#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class electrokineticCloud Declaration
\*---------------------------------------------------------------------------*/

class electrokineticCloud
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        electrokineticCloud(const electrokineticCloud&);

        //- Disallow default bitwise assignment
        void operator=(const electrokineticCloud&);


public:

    //- Runtime type information
    TypeName("electrokineticCloud");

    // Constructors

        //- Null constructor
        electrokineticCloud();


    //- Destructor
    virtual ~electrokineticCloud();


    // Member Functions

		//sources

			//- Return the particle volume fraction field
                //  Note: for particles belonging to this cloud only
                virtual const tmp<volScalarField> SPhi() const = 0;

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
