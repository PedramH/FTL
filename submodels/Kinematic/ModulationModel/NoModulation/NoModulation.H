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
    Foam::NoModulation

Description
    Place holder for 'none' option

\*---------------------------------------------------------------------------*/

#ifndef NoModulation_H
#define NoModulation_H

#include "ModulationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class NoModulation Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class NoModulation
:
    public ModulationModel<CloudType>
{
public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct from components
        NoModulation(const dictionary& dict, CloudType& owner);

        //- Construct copy
        NoModulation(NoModulation<CloudType>& dm);

        //- Construct and return a clone
        virtual autoPtr<ModulationModel<CloudType> > clone()
        {
            return autoPtr<ModulationModel<CloudType> >
            (
                new NoModulation<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~NoModulation();


    // Member Functions

        //- Update (disperse particles)
        virtual symmTensor update
        (
			const vector UcMean,
			const vector Uc,
			const vector Up0,
			const vector Upn,
			const vector Su
			//const typename CloudType::parcelType& p
		);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "NoModulation.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
