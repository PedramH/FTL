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
    Foam::MandoModulation

Description
    Class for particle turbulenc modulation model based on Mando's approach.

\*---------------------------------------------------------------------------*/

#ifndef MandoModulation_H
#define MandoModulation_H

#include "ModulationModel.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class MandoModulation Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class MandoModulation
:
    public ModulationModel<CloudType>
{
protected:

	// Protected data

		//- model coefficients
		Switch UpImplicit_;

public:

    //- Runtime type information
    TypeName("MandoModulation");


    // Constructors

        //- Construct from components
        MandoModulation(const dictionary& dict, CloudType& owner);

        //- Construct copy
        MandoModulation(MandoModulation<CloudType>& dm);

        //- Construct and return a clone
        virtual autoPtr<ModulationModel<CloudType> > clone()
        {
            return autoPtr<ModulationModel<CloudType> >
            (
                new MandoModulation<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~MandoModulation();


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
#   include "MandoModulation.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
