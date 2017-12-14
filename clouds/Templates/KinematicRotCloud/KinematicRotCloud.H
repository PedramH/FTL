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
    Foam::KinematicRotCloud

Description
    Templated base class for thermodynamic cloud

    - Adds to kinematic cloud
      - Particle rotation

SourceFiles
    KinematicRotCloudI.H
    KinematicRotCloud.C

\*---------------------------------------------------------------------------*/

#ifndef KinematicRotCloud_H
#define KinematicRotCloud_H

#include "KinematicCloud.H"
#include "kinematicRotCloud.H"

#include "ParticleTorqueList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

//template<class CloudType>
//class DeterminesticCollisionModel;

/*---------------------------------------------------------------------------*\
                        Class KinematicRotCloud Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class KinematicRotCloud
:
    public CloudType,
    public kinematicRotCloud
{
public:

    // Public typedefs

        //- Type of cloud this cloud was instantiated for
        typedef CloudType cloudType;

        //- Type of parcel the cloud was instantiated for
        typedef typename CloudType::particleType parcelType;

        //- Convenience typedef for this cloud type
        typedef KinematicRotCloud<CloudType> kinematicRotCloudType;

		//- Torque models type
        typedef ParticleTorqueList<KinematicRotCloud<CloudType> > torqueType;


private:

    // Private data

        //- Cloud copy pointer
        autoPtr<KinematicRotCloud<CloudType> > cloudCopyPtr_;


    // Private member functions

        //- Disallow default bitwise copy construct
        KinematicRotCloud(const KinematicRotCloud&);

        //- Disallow default bitwise assignment
        void operator=(const KinematicRotCloud&);


protected:

    // Protected data

        //- KinematicRot parcel constant properties
        typename parcelType::constantProperties constProps_;


        // References to the carrier gas fields

 			//- Curl of carrier phase velocity 
        	//autoPtr<volVectorField> curlUPtr_; //OmegaC memory


		//-Optional particle torques
        torqueType torques_;


        // References to the cloud sub-models

			//- statistics computation flags
			//bool stat2Comp_;


        // Reference to the particle integration schemes

            //- Angular Velocity integration
            autoPtr<vectorIntegrationScheme> OmegaIntegrator_;

		// Sources

            //- Momentum
            autoPtr<DimensionedField<vector, volMesh> > OmegaTrans_;

            //- Coefficient for carrier phase Omega equation
            //autoPtr<DimensionedField<scalar, volMesh> > OmegaCoeff_;


    // Protected Member Functions

         // Initialisation

            //- Set cloud sub-models
            void setModels();


        // Cloud evolution functions

            //- Reset state of cloud
            void cloudReset(KinematicRotCloud<CloudType>& c);


public:

    // Constructors

        //- Construct given carrier gas fields
        KinematicRotCloud
        (
			const word& cloudName,
            const volScalarField& rho,
            const volVectorField& U,
            const volScalarField& mu,
            const dimensionedVector& g,
            bool readFields = true
        );

        //- Copy constructor with new name
        KinematicRotCloud(KinematicRotCloud<CloudType>& c, const word& name);

        //- Copy constructor with new name - creates bare cloud
        KinematicRotCloud
        (
            const fvMesh& mesh,
            const word& name,
            const KinematicRotCloud<CloudType>& c
        );

        //- Construct and return clone based on (this) with new name
        virtual autoPtr<Cloud<parcelType> > clone(const word& name)
        {
            return autoPtr<Cloud<parcelType> >
            (
                new KinematicRotCloud(*this, name)
            );
        }

        //- Construct and return bare clone based on (this) with new name
        virtual autoPtr<Cloud<parcelType> > cloneBare(const word& name) const
        {
            return autoPtr<Cloud<parcelType> >
            (
                new KinematicRotCloud(this->mesh(), name, *this)
            );
        }


    //- Destructor
    virtual ~KinematicRotCloud();


    // Member Functions

		// Evaluation

            //- Cache fields
            //virtual void cacheFields(const bool store); //OmegaC memory

        // Access

            //- Return a reference to the cloud copy
            inline const KinematicRotCloud& cloudCopy() const;

            //- Return the constant properties
            inline const typename parcelType::constantProperties&
                constProps() const;

            //- Return access to the constant properties
            inline typename parcelType::constantProperties& constProps();

            //- Return const access to the carrier vorticity field
            //inline const volVectorField& curlU() const; //OmegaC memory


			//- Optional particle torques
            inline const torqueType& torques() const;

			//- Return the optional particle torques
            inline torqueType& torques();


            // Sub-models

 				//- Return the statistics computation flag
				//inline const bool& stat2Comp() const;

				//- Return the optional particle forces
				//inline bool& stat2Comp();


            // Integration schemes

                //-Return reference to angular velocity integration
                inline const vectorIntegrationScheme& OmegaIntegrator() const;

			 // Sources

                // Momentum

                    //- Return reference to momentum source
                    inline DimensionedField<vector, volMesh>& OmegaTrans();

                    //- Return const reference to momentum source
                    inline const DimensionedField<vector, volMesh>&
                        OmegaTrans() const;

                     //- Return coefficient for carrier phase U equation
                    //inline DimensionedField<scalar, volMesh>& OmegaCoeff();

                    //- Return const coefficient for carrier phase U equation
                    //inline const DimensionedField<scalar, volMesh>&
                    //    OmegaCoeff() const;


        // Check

            //- Maximum angular momentum
            inline scalar Omegamax() const;

            //- Minimum angular momentum
            inline scalar Omegamin() const;


        // Cloud evolution functions

            //- Set parcel thermo properties
            void setParcelThermoProperties
            (
                parcelType& parcel,
                const scalar lagrangianDt
            );

            //- Check parcel properties
            void checkParcelProperties
            (
                parcelType& parcel,
                const scalar lagrangianDt,
                const bool fullyDescribed
            );

            //- Store the current cloud state
            void storeState();

            //- Reset the current cloud to the previously stored state
            void restoreState();

            //- Reset the cloud source terms
            void resetSourceTerms();

            //- Apply relaxation to (steady state) cloud sources
            void relaxSources(const KinematicRotCloud<CloudType>& cloudOldTime);

            //- Apply scaling to (transient) cloud sources
            void scaleSources();

            //- Pre-evolve
            void preEvolve();

            //- Evolve the cloud
            void evolve();

			//- Particle motion
            template<class TrackData>
            void motion(TrackData& td);

			//- Particle field updates //added1
            //template<class TrackData> 
            //void preCollisionUpdates(TrackData& td);

			//- Post-evolve //change to public to be accessed by superclass kinematicCloud
            void postEvolve();


        // Mapping

            //- Remap the cells of particles corresponding to the
            //  mesh topology change with a default tracking data object
            virtual void autoMap(const mapPolyMesh&);


        // I-O

            //- Print cloud information
            void info();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "KinematicRotCloudI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "KinematicRotCloud.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
