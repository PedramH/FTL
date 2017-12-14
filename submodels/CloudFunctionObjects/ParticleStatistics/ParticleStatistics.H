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
    Foam::ParticleStatistics

Description
    Creates particle statistic to keep the total sum per cell and the average
    and rms dispersed phase velocity per cell. 

    Add to particlePropertiesDict 

    cloudFunctions
    {
      particleStatistics1
      {
		 type			 particleStatistics;
         resetFields      false;          // on a restart, reset the last fields or continue

		 gNum		     2;               // Number of size groups
		 dgMin			 (4e-6 7e-6);     // the minimum diameter of each group
		 dgMax			 (7e-6 15e-6);    // the maximum diameter of each group

		 UMean			  true;			  // calculate mean velocity
         uPrime2Mean      false;          // Also keep and write the ux'2, uy'2 and uz'2, uxuy, uxuz, uyuz
         d10			  false;          // mean diameter
		 d32			  false;          // Sauter diameter
         timeStart        0.0;            // Start averaging time
         timeEnd          1.0;            // Stop averaging time
      }
   }

Eelco van Vliet
Tata Steel, Sep 2012 :Original
Ehsan Amani, June 2015 :Extended


SourceFiles
    ParticleStatistics.C

\*---------------------------------------------------------------------------*/

#ifndef ParticleStatistics_H
#define ParticleStatistics_H

#include "CloudFunctionObject.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class ParticleStatistics Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class ParticleStatistics
:
    public CloudFunctionObject<CloudType>
{
    // Private Data

        // Typedefs

            //- Convenience typedef for parcel type
            typedef typename CloudType::parcelType parcelType;


        //- Void fraction field
        List<autoPtr<volScalarField> > ParticleSumPtr_;
		List<autoPtr<volScalarField> > ParticleCSumPtr_;
        List<autoPtr<volVectorField> > ParticleUSumPtr_;
        List<autoPtr<volSymmTensorField> > ParticleUsqSumPtr_;
		List<autoPtr<volScalarField> > Particled10SumPtr_;
		List<autoPtr<volScalarField> > Particled2SumPtr_;
		List<autoPtr<volScalarField> > Particled3SumPtr_;

        // initialise fields to zero on start averaging
        Switch resetFields_; 

		// group number
        label gNum_;

		// calculate mean diameter statistics
        Switch UMean_; 

        // calculate uprime statistics
        Switch uPrime2Mean_; 

        // write the particle PDF
        Switch writed10_; 

		// write the particle PDF
        Switch writed32_; 

        // start averaging at
        scalar timeStart_; 

        // end averaging at
        scalar timeEnd_; 

		// the minimum diameter of each group
        List<scalar> dgMin_;

		// the maximum diameter of each group
		List<scalar> dgMax_;	

		virtual void nullifyPtrLists();

protected:

    // Protected Member Functions

        //- Write post-processing info
        virtual void write();


public:

    //- Runtime type information
    TypeName("particleStatistics");


    // Constructors

        //- Construct from dictionary
        ParticleStatistics
		(
            const dictionary& dict,
            CloudType& owner,
            const word& modelName
        );

        //- Construct copy
        ParticleStatistics(const ParticleStatistics<CloudType>& ps);

        //- Construct and return a clone
        virtual autoPtr<CloudFunctionObject<CloudType> > clone() const
        {
            return autoPtr<CloudFunctionObject<CloudType> >
            (
                new ParticleStatistics<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~ParticleStatistics();


    // Member Functions
    //
    // Access
      
        //- Return const access to the resetFields flag
        inline const Switch& resetFields() const;

		//- Return const access to the resetFields flag
        inline const int& gNum() const;

		//- Return const access to the resetFields flag
        inline const Switch& UMean() const;

        //- Return const access to the resetFields flag
        inline const Switch& uPrime2Mean() const;

        //- Return const access to the resetFields flag
        inline const Switch& writed10() const;

		//- Return const access to the resetFields flag
        inline const Switch& writed32() const;

        //- Return start avering time 
        inline label timeStart() const;

        //- Return end avering time 
        inline label timeEnd() const;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

        // Evaluation

            //- Pre-evolve hook
            virtual void preEvolve();

            //- Post-evolve hook
            virtual void postEvolve();

            //- Post-move hook
            virtual void postMove
            (   
                const parcelType& p,
                const label cellI,
                const scalar dt,
				const point& position0,
				bool& keepParticle
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//
#include "ParticleStatisticsI.H"

#ifdef NoRepository
#   include "ParticleStatistics.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
