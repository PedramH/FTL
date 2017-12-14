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
    Foam::mykinematicParcelInjectionData

Description
    Container class to provide injection data for mykinematicdynamic parcels

SourceFiles
    mykinematicParcelInjectionData.C

\*---------------------------------------------------------------------------*/

#ifndef mykinematicParcelInjectionData_H
#define mykinematicParcelInjectionData_H

#include "kinematicParcelInjectionData.H"
#include "symmTensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class mykinematicParcelInjectionData;

// Forward declaration of friend functions

Ostream& operator<<
(
    Ostream&,
    const mykinematicParcelInjectionData&
);

Istream& operator>>
(
    Istream&,
    mykinematicParcelInjectionData&
);

/*---------------------------------------------------------------------------*\
                 Class mykinematicParcelInjectionData Declaration
\*---------------------------------------------------------------------------*/

class mykinematicParcelInjectionData
:
    public kinematicParcelInjectionData
{
protected:

    // Parcel properties

		 //- diameter interval length [m]
        scalar deltaD_;

        //- (particle) Reynolds stress [m2/s2]
        symmTensor uiuj_;

		//- fluid phase kinetic energy at injection location [m2/s2]
		scalar kf_;


public:

    //- Runtime type information
    TypeName("mykinematicParcelInjectionData");

    // Constructors

        //- Null constructor
        mykinematicParcelInjectionData();

        //- Construct from dictionary
        mykinematicParcelInjectionData(const dictionary& dict);

        //- Construct from Istream
        mykinematicParcelInjectionData(Istream& is);


    //-Destructor
    virtual ~mykinematicParcelInjectionData();


    // Access

		//- Return const access to deltaD
        inline scalar deltaD() const;

        //- Return const access to the RS
        inline symmTensor uiuj() const;

		//- Return const access to the kf
        inline scalar kf() const;


    // Edit

		//- Return const access to deltaD
        inline scalar& deltaD();

        //- Return access to the RS
        inline symmTensor& uiuj();

		//- Return access to the kf
        inline scalar& kf();


    // I-O

        // Friend Operators

            friend bool operator==
            (
                const mykinematicParcelInjectionData& a,
                const mykinematicParcelInjectionData& b
            )
            {
                notImplemented
                (
                    "operator=="
                    "("
                        "const mykinematicParcelInjectionData&, "
                        "const mykinematicParcelInjectionData&"
                    ")"
                );

                return false;
            }

            friend bool operator!=
            (
                const mykinematicParcelInjectionData& a,
                const mykinematicParcelInjectionData& b
            )
            {
                notImplemented
                (
                    "operator=="
                    "("
                        "const mykinematicParcelInjectionData&, "
                        "const mykinematicParcelInjectionData&"
                    ")"
                );

                return false;
            }

            //- Ostream operator
            friend Ostream& operator<<
            (
                Ostream& os,
                const mykinematicParcelInjectionData& data
            );

            //- Istream operator
            friend Istream& operator>>
            (
                Istream& is,
                mykinematicParcelInjectionData& data
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "mykinematicParcelInjectionDataI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
