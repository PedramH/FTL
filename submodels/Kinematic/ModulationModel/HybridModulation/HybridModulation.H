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
    Foam::HybridModulation

Description
    Class for particle turbulenc modulation model based on hybrid approach.

\*---------------------------------------------------------------------------*/

#ifndef HybridModulation_H
#define HybridModulation_H

#include "ModulationModel.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class HybridModulation Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class HybridModulation
:
    public ModulationModel<CloudType>
{
protected:

	// Protected data

		//- model coefficients
		scalar alpha_;

		scalar beta_;

		Switch UpImplicit_;


public:

    //- Runtime type information
    TypeName("hybridModulation");


    // Constructors

        //- Construct from components
        HybridModulation(const dictionary& dict, CloudType& owner);

        //- Construct copy
        HybridModulation(HybridModulation<CloudType>& dm);

        //- Construct and return a clone
        virtual autoPtr<ModulationModel<CloudType> > clone()
        {
            return autoPtr<ModulationModel<CloudType> >
            (
                new HybridModulation<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~HybridModulation();


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
#   include "HybridModulation.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
