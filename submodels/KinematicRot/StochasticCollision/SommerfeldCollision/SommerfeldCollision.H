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
    Foam::SommerfeldCollision

Description
    Collision model by 


\*---------------------------------------------------------------------------*/

#ifndef SommerfeldCollision_H
#define SommerfeldCollision_H

#include "StochasticCollisionModel.H"
//#include "liquidMixtureProperties.H"

#include "Random.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                      Class SommerfeldCollision Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class SommerfeldCollision
:
    public StochasticCollisionModel<CloudType>
{
protected:

    // Protected Data

        //- Convenience typedef to the cloud's parcel type
        typedef typename CloudType::parcelType parcelType;

        //const liquidMixtureProperties& liquids_;

        //- Coalescence activation switch
        //Switch coalescence_;

		//- Random number generator
		Random ranGen_;

		//- restituation, e
        scalar e_;

        //- friction factor, mu
        scalar mu_;

		//activating repulsive correction
		Switch repulsiveCorrection_;

		//- Particle mean Velocity
		const DimensionedField<vector, volMesh>* UpMean_;

		//- Particle Velocity varience
		const DimensionedField<vector, volMesh>* upup_;

		//- Particle mean number density
		const DimensionedField<scalar, volMesh>* NpMean_;

		//- Particle mean diameter
		const DimensionedField<scalar, volMesh>* dpMean_;

		//- Particle diameter varience
		const DimensionedField<scalar, volMesh>* dpdp_;

		//- Particle mean angular Velocity
		const DimensionedField<vector, volMesh>* OmegapMean_;

		// Locally cached carrier fields

            //- Turbulence k
            const volScalarField* kPtr_;

            //- Take ownership of the k field
            bool ownK_;

            //- Turbulence epsilon
            const volScalarField* epsilonPtr_;

            //- Take ownership of the epsilon field
            bool ownEpsilon_;

			//- rho
            const volScalarField* rhoPtr_;

            //- Take ownership of the rho field
            bool ownRho_;

			//- rho
            const volScalarField* muPtr_;

            //- Take ownership of the mu field
            bool ownMu_;


    // Protected Member Functions

        //- Main collision routine
        virtual void collide(const scalar dt);

        //- Collide parcels and return true if mass has changed
        virtual bool collideParcels
        (
            const scalar dt,
            parcelType& p
        );

        // L is the larger drop and S is the smaller
        virtual scalar limitingDistance 
        (
            const scalar muc,
			const scalar rhoc,
			const scalar dL,
			const vector UL,
			const scalar rhoS,
			const scalar dS,
			const vector US
        );

		//- Return the k field from the turbulence model
        tmp<volScalarField> kModel() const;

        //- Return the epsilon field from the turbulence model
        tmp<volScalarField> epsilonModel() const;

		//- Return the rho field from the continous phase
        tmp<volScalarField> rhoModel() const;

        //- Return the mu field from the continous phase
        tmp<volScalarField> muModel() const;


public:

    //- Runtime type information
    TypeName("SommerfeldCollision");


    // Constructors

        //- Construct from dictionary
        SommerfeldCollision
        (
            const dictionary& dict,
            CloudType& cloud,
            const word& modelName = typeName
        );

        //- Construct copy
        SommerfeldCollision(SommerfeldCollision<CloudType>& cm);

        //- Construct and return a clone
        virtual autoPtr<StochasticCollisionModel<CloudType> > clone() //const
        {
            return autoPtr<StochasticCollisionModel<CloudType> >
            (
                new SommerfeldCollision<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~SommerfeldCollision();


	// Member Functions

		//- Cache carrier fields
        virtual void cacheFields(const bool store);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "SommerfeldCollision.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
