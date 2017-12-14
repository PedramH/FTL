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
    Foam::KinematicRotParcel

Description
    Thermodynamic parcel class with one/two-way coupling with the continuous
    phase. Includes Kinematic parcel sub-models, plus:
    - heat transfer

SourceFiles
    KinematicRotParcelI.H
    KinematicRotParcel.C
    KinematicRotParcelIO.C

\*---------------------------------------------------------------------------*/

#ifndef KinematicRotParcel_H
#define KinematicRotParcel_H

#include "particle.H"
#include "demandDrivenEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class ParcelType>
class KinematicRotParcel;

template<class ParcelType>
Ostream& operator<<
(
    Ostream&,
    const KinematicRotParcel<ParcelType>&
);

/*---------------------------------------------------------------------------*\
                       Class KinematicRotParcel Declaration
\*---------------------------------------------------------------------------*/

template<class ParcelType>
class KinematicRotParcel
:
    public ParcelType
{
public:

    //- Class to hold kinematicRot particle constant properties
    class constantProperties
    :
        public ParcelType::constantProperties
    {

        // Private data

            //- Particle initial angular velocity [rad/s]
            demandDrivenEntry<vector> Omega0_;


    public:

        // Constructors

            //- Null constructor
            constantProperties();

            //- Copy constructor
            constantProperties(const constantProperties& cp);

            //- Construct from dictionary
            constantProperties(const dictionary& parentDict);


        // Member functions

            // Access

                //- Return const access to the particle initial angular velocity [rad/s]
                inline vector Omega0() const;
    };


    template<class CloudType>
    class TrackingData
    :
        public ParcelType::template TrackingData<CloudType>
    {
    private:

        // Private data

				// Particle Statistics

				//- CRW dispersion model data calculation
		        //Switch particleStatistics1_;

				//- Sommerfeld's collision model data calculation
		        //Switch particleStatistics2_; //added1

				//- Velocity average
                //autoPtr<DimensionedField<vector, volMesh> > UPStat_;

				//- Velocity component squared average
                autoPtr<DimensionedField<vector, volMesh> > OmegaPStat_;


            // Interpolators for continuous phase fields

                //- CurlU field interpolator
                //autoPtr<interpolation<vector> > curlUInterp_; //OmegaC memory

    public:

        typedef typename ParcelType::template TrackingData<CloudType>::trackPart
            trackPart;
		typedef typename ParcelType::template TrackingData<CloudType>::fieldUpdate
            fieldUpdate; //added1

        // Constructors

            //- Construct from components
            inline TrackingData
            (
                CloudType& cloud,
                trackPart part = ParcelType::template
                    TrackingData<CloudType>::tpLinearTrack,
				fieldUpdate updates = ParcelType::template //added1
                    TrackingData<CloudType>::motion
            );

		//- Destructor
    	virtual ~TrackingData(); //added1

		
		//- Update the particle averages
        virtual inline void updateParticleStatistics(); //(CloudType& cloud); //added1


        // Member functions

            //- Return const access to the interpolator for continuous
            //  phase vorticity field
            //inline const interpolation<vector>& curlUInterp() const; //OmegaC memory
    };


protected:

    // Protected data

        // Parcel properties

            //Angular Velocity of Parcel [rad/s]
            vector Omega_;

		// Cell-based quantities

            //- Angular velocity [rad/s]
            vector Omegac_;


    // Protected Member Functions

        //- Calculate new particle temperature
        template<class TrackData>
        const vector calcOmega
        (
            TrackData& td,
            const scalar dt,           // timestep
            const label cellI,         // owner cell
            const scalar Re,           // Reynolds number
            const scalar mu,           // local carrier viscosity
            const scalar momentOfInertia,         // mass
            const vector& Somega,          // explicit particle angular momentum source
            vector& dOmegaTrans,           // angular momentum transfer to carrier
            scalar& Spomega                // linearised drag coefficient
        ); //const;


public:

    // Static data members

        //- Runtime type information
        TypeName("KinematicRotParcel");

        //- String representation of properties
        AddToPropertyList
        (
            ParcelType,
            " (Omegax Omegay Omegaz)"
        );


    // Constructors

        //- Construct from owner, position, and cloud owner
        //  Other properties initialised as null
        inline KinematicRotParcel
        (
            const polyMesh& mesh,
            const vector& position,
            const label cellI,
            const label tetFaceI,
            const label tetPtI
        );

        //- Construct from components
        inline KinematicRotParcel
        (
            const polyMesh& mesh,
            const vector& position,
            const label cellI,
            const label tetFaceI,
            const label tetPtI,
            const label typeId,
            const scalar nParticle0,
            const scalar d0,
            const scalar dTarget0,
            const vector& U0,
            const constantProperties& constProps
        );

        //- Construct from Istream
        KinematicRotParcel
        (
            const polyMesh& mesh,
            Istream& is,
            bool readFields = true
        );

        //- Construct as a copy
        KinematicRotParcel(const KinematicRotParcel& p);

        //- Construct as a copy
        KinematicRotParcel(const KinematicRotParcel& p, const polyMesh& mesh);

        //- Construct and return a (basic particle) clone
        virtual autoPtr<particle> clone() const
        {
            return autoPtr<particle>(new KinematicRotParcel(*this));
        }

        //- Construct and return a (basic particle) clone
        virtual autoPtr<particle> clone(const polyMesh& mesh) const
        {
            return autoPtr<particle>(new KinematicRotParcel(*this, mesh));
        }

        //- Factory class to read-construct particles used for
        //  parallel transfer
        class iNew
        {
            const polyMesh& mesh_;

        public:

            iNew(const polyMesh& mesh)
            :
                mesh_(mesh)
            {}

            autoPtr<KinematicRotParcel<ParcelType> > operator()(Istream& is) const
            {
                return autoPtr<KinematicRotParcel<ParcelType> >
                (
                    new KinematicRotParcel<ParcelType>(mesh_, is, true)
                );
            }
        };


    // Member Functions

        // Access

            //- Return access to Angular velocity
            inline const vector& Omega() const; 

			//- Return const access to carrier angular velocity
            inline vector Omegac() const;


        // Edit

            //- Return access to Angular velocity
            inline vector& Omega(); 


		// Helper functions

            //- Cell owner moment of Inerta 
            inline scalar momentOfInertiaCell(const label cellI) const;


        // Main calculation loop

            //- Set cell values
            template<class TrackData>
            void setCellValues
            (
                TrackData& td,
                const scalar dt,
                const label cellI
            );

            //- Correct cell values using latest transfer information
            template<class TrackData>
            void cellValueSourceCorrection
            (
                TrackData& td,
                const scalar dt,
                const label cellI
            );

            //- Update parcel properties over the time interval
            template<class TrackData>
            void calc
            (
                TrackData& td,
                const scalar dt,
                const label cellI
            );

		
		// Patch interactions

			//- Transform the physical properties of the particle
            //  according to the given transformation tensor
            virtual void transformProperties(const tensor& T);

            //- Transform the physical properties of the particle
            //  according to the given separation vector
            virtual void transformProperties(const vector& separation);


        // I-O

            //- Read
            template<class CloudType>
            static void readFields(CloudType& c);

            //- Write
            template<class CloudType>
            static void writeFields(const CloudType& c);


    // Ostream Operator

        friend Ostream& operator<< <ParcelType>
        (
            Ostream&,
            const KinematicRotParcel<ParcelType>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "KinematicRotParcelI.H"
#include "KinematicRotParcelTrackingDataI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "KinematicRotParcel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
