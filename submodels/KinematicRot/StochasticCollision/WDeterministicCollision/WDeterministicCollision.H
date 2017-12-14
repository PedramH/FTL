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
    Foam::WDeterministicCollision

Description
    Collision model by 


\*---------------------------------------------------------------------------*/

#ifndef WDeterministicCollision_H
#define WDeterministicCollision_H

#include "StochasticCollisionModel.H"
//#include "liquidMixtureProperties.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                      Class WDeterministicCollision Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class WDeterministicCollision
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

		//- restituation, e
        scalar e_;

        //- friction factor, mu
        scalar mu_;

		//activating repulsive correction
		Switch repulsiveCorrection_;

	// Locally cached carrier fields

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
            parcelType& p1,
            parcelType& p2,
            scalar& m1,
            scalar& m2
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

		//- Return the rho field from the continous phase
        tmp<volScalarField> rhoModel() const;

        //- Return the mu field from the continous phase
        tmp<volScalarField> muModel() const;


public:

    //- Runtime type information
    TypeName("WDeterministicCollision");


    // Constructors

        //- Construct from dictionary
        WDeterministicCollision
        (
            const dictionary& dict,
            CloudType& cloud,
            const word& modelName = typeName
        );

        //- Construct copy
        WDeterministicCollision(WDeterministicCollision<CloudType>& cm);

        //- Construct and return a clone
        virtual autoPtr<StochasticCollisionModel<CloudType> > clone() //const
        {
            return autoPtr<StochasticCollisionModel<CloudType> >
            (
                new WDeterministicCollision<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~WDeterministicCollision();


	// Member Functions

		//- Cache carrier fields
        virtual void cacheFields(const bool store);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "WDeterministicCollision.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
