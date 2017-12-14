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
    Foam::ElectrokineticCloud

Description
    Templated base class for thermodynamic cloud

    - Adds to kinematic cloud
      - Particle rotation

SourceFiles
    ElectrokineticCloudI.H
    ElectrokineticCloud.C

\*---------------------------------------------------------------------------*/

#ifndef ElectrokineticCloud_H
#define ElectrokineticCloud_H

#include "electrokineticCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

//template<class CloudType>
//class DeterminesticCollisionModel;

/*---------------------------------------------------------------------------*\
                        Class ElectrokineticCloud Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class ElectrokineticCloud
:
    public CloudType,
    public electrokineticCloud
{
public:

    // Public typedefs

        //- Type of cloud this cloud was instantiated for
        typedef CloudType cloudType;

        //- Type of parcel the cloud was instantiated for
        typedef typename CloudType::particleType parcelType;

        //- Convenience typedef for this cloud type
        typedef ElectrokineticCloud<CloudType> electrokineticCloudType;


private:

    // Private data

        //- Cloud copy pointer
        autoPtr<ElectrokineticCloud<CloudType> > cloudCopyPtr_;


    // Private member functions

        //- Disallow default bitwise copy construct
        ElectrokineticCloud(const ElectrokineticCloud&);

        //- Disallow default bitwise assignment
        void operator=(const ElectrokineticCloud&);


protected:

    // Protected data

        //- Electrokinetic parcel constant properties
        typename parcelType::constantProperties constProps_;


        // References to the carrier gas fields

 			//- Velocity [m/s]
            const volVectorField& E_;


    // Protected Member Functions

         // Initialisation

            //- Set cloud sub-models
            void setModels();


        // Cloud evolution functions

            //- Reset state of cloud
            void cloudReset(ElectrokineticCloud<CloudType>& c);


public:

    // Constructors

        //- Construct given carrier gas fields
        ElectrokineticCloud
        (
			const word& cloudName,
            const volScalarField& rho,
            const volVectorField& U,
            const volScalarField& mu,
			const volVectorField& E,
            const dimensionedVector& g,
            bool readFields = true
        );

        //- Copy constructor with new name
        ElectrokineticCloud(ElectrokineticCloud<CloudType>& c, const word& name);

        //- Copy constructor with new name - creates bare cloud
        ElectrokineticCloud
        (
            const fvMesh& mesh,
            const word& name,
            const ElectrokineticCloud<CloudType>& c
        );

        //- Construct and return clone based on (this) with new name
        virtual autoPtr<Cloud<parcelType> > clone(const word& name)
        {
            return autoPtr<Cloud<parcelType> >
            (
                new ElectrokineticCloud(*this, name)
            );
        }

        //- Construct and return bare clone based on (this) with new name
        virtual autoPtr<Cloud<parcelType> > cloneBare(const word& name) const
        {
            return autoPtr<Cloud<parcelType> >
            (
                new ElectrokineticCloud(this->mesh(), name, *this)
            );
        }


    //- Destructor
    virtual ~ElectrokineticCloud();


    // Member Functions

		// Evaluation

            //- Cache fields
            //virtual void cacheFields(const bool store); //OmegaC memory

        // Access

            //- Return a reference to the cloud copy
            inline const ElectrokineticCloud& cloudCopy() const;

            //- Return the constant properties
            inline const typename parcelType::constantProperties&
                constProps() const;

            //- Return access to the constant properties
            inline typename parcelType::constantProperties& constProps();

            //- Return const access to the carrier electric field
            inline const volVectorField& E() const; 


		// Fields

                //- Return the potential field source term
                //  Note: for particles belonging to this cloud only
                inline const tmp<volScalarField> SPhi() const;


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
            void relaxSources(const ElectrokineticCloud<CloudType>& cloudOldTime);

            //- Apply scaling to (transient) cloud sources
            void scaleSources();

            //- Pre-evolve
            void preEvolve();

            //- Evolve the cloud
            void evolve();

			//- Particle motion
            template<class TrackData>
            void motion(TrackData& td);

			//- Particle field updates before stochastic collision //added1
            /*template<class TrackData>
            void preCollisionUpdates(TrackData& td);*/

			//- Post-evolve //change to public to be accessed by superclass electrokineticCloud
            void postEvolve();


        // Mapping

            //- Remap the cells of particles corresponding to the
            //  mesh topology change with a default tracking data object
            virtual void autoMap(const mapPolyMesh&);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ElectrokineticCloudI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "ElectrokineticCloud.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
