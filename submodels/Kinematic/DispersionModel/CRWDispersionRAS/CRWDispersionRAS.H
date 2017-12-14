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
    Foam::CRWDispersionRAS

Description
    The velocity is perturbed in random direction, with a
    Gaussian random number distribution with variance sigma.
    where sigma is defined below

\*---------------------------------------------------------------------------*/

#ifndef CRWDispersionRAS_H
#define CRWDispersionRAS_H

#include "DispersionRASModel.H"
#include "Random.H"
//#include "wallDist.H"
#include "wallDistReflection.H"
#include "interpolation.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class CRWDispersionRAS Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class CRWDispersionRAS
:
    public DispersionRASModel<CloudType>
{
protected:

	// Protected data

		//- Random number generator
		Random ranGen_;

		//- model coefficients
		scalar CL_;

		scalar beta_;

		Switch anisotropy_;

		scalar uiuiMin_;

		Switch nearWallT_;

		scalar Cmu_;

		scalar yPlusBL_;

		//- Rf 
            const volSymmTensorField* RfPtr_;

            //- Take ownership of the Rf
            bool ownRf_;

		//- drift 
            const volVectorField* driftPtr_;

            //- Take ownership of the drift
            bool ownDrift_;

		//- Particle mean Velocity
		const DimensionedField<vector, volMesh>* UpMean_;

		//wallDist y_;
		wallDistReflection yr_;


		// Protected Functions

			//- Return the Rf field from the turbulence model
		    tmp<volSymmTensorField> RfModel();

		    //- Return the drift field from the turbulence model
		    tmp<volVectorField> driftModel();


private:

		//- Density interpolator
                autoPtr<interpolation<scalar> > kInterp_;

                //- Velocity interpolator
                autoPtr<interpolation<scalar> > epsInterp_;

				//- Dynamic viscosity interpolator
                autoPtr<interpolation<symmTensor> > RfInterp_;

                //- Dynamic viscosity interpolator
                autoPtr<interpolation<vector> > driftInterp_;

		//- wall normal vector interpolator
        autoPtr<interpolation<vector> > nwInterp_;

		//- wall distance interpolator
        autoPtr<interpolation<scalar> > yInterp_;

	
public:

    //- Runtime type information
    TypeName("CRWDispersionRAS");


    // Constructors

        //- Construct from components
        CRWDispersionRAS(const dictionary& dict, CloudType& owner);

        //- Construct copy
        CRWDispersionRAS(CRWDispersionRAS<CloudType>& dm);

        //- Construct and return a clone
        virtual autoPtr<DispersionModel<CloudType> > clone()
        {
            return autoPtr<DispersionModel<CloudType> >
            (
                new CRWDispersionRAS<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~CRWDispersionRAS();


    // Member Functions

		//- Cache carrier fields
        virtual void cacheFields(const bool store);

		// I-O

            //- Write
            virtual void write(Ostream& os) const;

        //- Update (disperse particles)
        virtual vector update
        /*(
            const scalar dt,
            const label cellI,
            const vector& U,
            const vector& Uc,
            vector& UTurb,
            scalar& tTurb
        );*/
        (
            const scalar dt,
            const label cellI,
            const vector& U,
            const vector& Uc,
            vector& UTurb,
            scalar& tTurb,
			const scalar& muc,
			const scalar& rhoc,
			const typename CloudType::parcelType& p
        );

		virtual void tTurbUpdate
        (
			scalar& tTurb,
			scalar& taop
		);

		void CRW_Dispersion
		(
			const scalar dt,
			const label cellI,
			const vector& U,
			const vector& Uc,
			vector& UTurb,
			scalar& taop,
			const scalar& k,
			const scalar& epsilon,
			symmTensor Rf,
			vector driftCorr,
			const vector UpMean
		);

		tensor CholskyDecompose(symmTensor& sisj);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "CRWDispersionRAS.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
