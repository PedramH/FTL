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
    Foam::ElectrokineticParcel

Description
    Thermodynamic parcel class with one/two-way coupling with the continuous
    phase. Includes Kinematic parcel sub-models, plus:
    - electrokinetic

SourceFiles
    ElectrokineticParcelI.H
    ElectrokineticParcel.C
    ElectrokineticParcelIO.C

\*---------------------------------------------------------------------------*/

#ifndef ElectrokineticParcel_H
#define ElectrokineticParcel_H

#include "particle.H"
#include "demandDrivenEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class ParcelType>
class ElectrokineticParcel;

template<class ParcelType>
Ostream& operator<<
(
    Ostream&,
    const ElectrokineticParcel<ParcelType>&
);

/*---------------------------------------------------------------------------*\
                       Class ElectrokineticParcel Declaration
\*---------------------------------------------------------------------------*/

template<class ParcelType>
class ElectrokineticParcel
:
    public ParcelType
{
public:

	//- Class to hold electrokinetic particle constant properties
    class constantProperties
    :
        public ParcelType::constantProperties
    {

        // Private data

            //- Particle initial valance number [-]
            demandDrivenEntry<scalar> z0_;


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

                //- Return const access to the particle valance number [-]
                inline scalar z0() const;
    };


    template<class CloudType>
    class TrackingData
    :
        public ParcelType::template TrackingData<CloudType>
    {
    private:

        // Private data

            // Interpolators for continuous phase fields

                //- Electric field interpolator
                autoPtr<interpolation<vector> > EInterp_; //electrokineticFoam

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

            //- Return conat access to the interpolator for continuous
            //  phase electric field
            inline const interpolation<vector>& EInterp() const; 
    };


protected:

    // Protected data

        // Parcel properties

            //- valance number of ion [s]
            scalar z_; 

		// Cell-based quantities

            //- Electric field [m/s]
            vector Ec_; 


    // Protected Member Functions


public:

    // Static data members

        //- Runtime type information
        TypeName("ElectrokineticParcel");

        //- String representation of properties
        AddToPropertyList
        (
            ParcelType,
            " valance" 
        );


    // Constructors

        //- Construct from owner, position, and cloud owner
        //  Other properties initialised as null
        inline ElectrokineticParcel
        (
            const polyMesh& mesh,
            const vector& position,
            const label cellI,
            const label tetFaceI,
            const label tetPtI
        );

        //- Construct from components
        inline ElectrokineticParcel
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
        ElectrokineticParcel
        (
            const polyMesh& mesh,
            Istream& is,
            bool readFields = true
        );

        //- Construct as a copy
        ElectrokineticParcel(const ElectrokineticParcel& p);

        //- Construct as a copy
        ElectrokineticParcel(const ElectrokineticParcel& p, const polyMesh& mesh);

        //- Construct and return a (basic particle) clone
        virtual autoPtr<particle> clone() const
        {
            return autoPtr<particle>(new ElectrokineticParcel(*this));
        }

        //- Construct and return a (basic particle) clone
        virtual autoPtr<particle> clone(const polyMesh& mesh) const
        {
            return autoPtr<particle>(new ElectrokineticParcel(*this, mesh));
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

            autoPtr<ElectrokineticParcel<ParcelType> > operator()(Istream& is) const
            {
                return autoPtr<ElectrokineticParcel<ParcelType> >
                (
                    new ElectrokineticParcel<ParcelType>(mesh_, is, true)
                );
            }
        };


    // Member Functions

        // Access

            //- Return access to valance number
            inline const scalar& z() const; 

			//- Return const access to carrier electric field
            inline vector Ec() const;


        // Edit

            //- Return access to valance number
            inline scalar& z(); 


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
            const ElectrokineticParcel<ParcelType>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ElectrokineticParcelI.H"
#include "ElectrokineticParcelTrackingDataI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ElectrokineticParcel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
