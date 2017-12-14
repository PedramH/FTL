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
    Foam::frontTrackingCloud

Description
    Virtual abstract base class for templated frontTrackingCloud

SourceFiles
    frontTrackingCloud.H

\*---------------------------------------------------------------------------*/

#ifndef frontTrackingCloud_H
#define frontTrackingCloud_H

#include "typeInfo.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class frontTrackingCloud Declaration
\*---------------------------------------------------------------------------*/

class frontTrackingCloud
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        frontTrackingCloud(const frontTrackingCloud&);

        //- Disallow default bitwise assignment
        void operator=(const frontTrackingCloud&);


public:

    //- Runtime type information
    TypeName("frontTrackingCloud");

    // Constructors

        //- Null constructor
        frontTrackingCloud();


    // Member functions

        // Check

            //-  Number of parcels
            virtual label nParcels() const = 0;

            //- Total mass in system
            virtual scalar massInSystem() const = 0;

            //- Total linear momentum of the system
            //virtual vector linearMomentumOfSystem() const = 0;

            //- Total linear kinetic energy in the system
            //virtual scalar linearKineticEnergyOfSystem() const = 0;

            //- Total rotational kinetic energy in the system
//            virtual scalar rotationalKineticEnergyOfSystem() const = 0;

            //- Penetration for percentage of the current total mass
//            virtual scalar penetration(const scalar& fraction) const = 0;

            //- Mean diameter Dij
            //virtual scalar Dij(const label i, const label j) const = 0;

            //- Max diameter
            //virtual scalar Dmax() const = 0;


            // Fields

                //- Volume swept rate of parcels per cell
                virtual const tmp<volScalarField> vDotSweep() const = 0;

                //- Return the particle volume fraction field
                //  Note: for particles belonging to this cloud only
                virtual const tmp<volScalarField> theta() const = 0;

                //- Return the particle mass fraction field
                //  Note: for particles belonging to this cloud only
                virtual const tmp<volScalarField> alpha() const = 0;

                //- Return the particle effective density field
                //  Note: for particles belonging to this cloud only
                virtual const tmp<volScalarField> rhoEff() const = 0;


    //- Destructor
    virtual ~frontTrackingCloud();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
