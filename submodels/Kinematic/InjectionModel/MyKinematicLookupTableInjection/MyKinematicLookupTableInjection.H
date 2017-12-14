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
    Foam::MyKinematicLookupTableInjection

Description
    Particle injection sources read from look-up table. Each row corresponds to
    an injection site.

    (
        (x y z) (u v w) d rho mDot deltaD (uu uv uw vv vw ww) kf  // injector 1
        (x y z) (u v w) d rho mDot deltaD (uu uv uw vv vw ww) kf  // injector 2
        ...
        (x y z) (u v w) d rho mDot deltaD (uu uv uw vv vw ww) kf  // injector N
    );

    where:
        x, y, z = global cartesian co-ordinates [m]
        u, v, w = global cartesian velocity components [m/s]
		d       = (mean) diameter [m]
        rho     = density [kg/m3]
        mDot    = mass flow rate [kg/m3]
		deltaD  = diameter range (d-deltaD/2,d+deltaD/2) (uniform distribution) [m]
		uu ...  = global Reynolds stress components [m/s]
		kf      = fluid phase kinetic energy at injection location [m2/s2]
		

SourceFiles
    MyKinematicLookupTableInjection.C

\*---------------------------------------------------------------------------*/

#ifndef MyKinematicLookupTableInjection_H
#define MyKinematicLookupTableInjection_H

#include "InjectionModel.H"
#include "mykinematicParcelInjectionDataIOList.H"
#include "Random.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class MyKinematicLookupTableInjection Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class MyKinematicLookupTableInjection
:
    public InjectionModel<CloudType>
{
    // Private data

		//- Random number generator
		Random ranGen_;

        //- Name of file containing injector/parcel data
        const word inputFileName_;

        //- Injection duration - common to all injection sources
        scalar duration_;

        //- Number of parcels per injector - common to all injection sources
        const scalar parcelsPerSecond_;

        //- Flag to indicate to randomise injection positions
        bool randomise_;

        //- List of injectors
        mykinematicParcelInjectionDataIOList injectors_;

        //- List of cell labels corresoponding to injector positions
        labelList injectorCells_;

        //- List of tetFace labels corresoponding to injector positions
        labelList injectorTetFaces_;

        //- List of tetPt labels corresoponding to injector positions
        labelList injectorTetPts_;


public:

    //- Runtime type information
    TypeName("myKinematicLookupTableInjection");


    // Constructors

        //- Construct from dictionary
        MyKinematicLookupTableInjection
        (
            const dictionary& dict,
            CloudType& owner,
            const word& modelName
        );

        //- Construct copy from owner cloud and injection model
        MyKinematicLookupTableInjection
        (
            const MyKinematicLookupTableInjection<CloudType>& im
        );

        //- Construct and return a clone using supplied owner cloud
        virtual autoPtr<InjectionModel<CloudType> > clone() const
        {
            return autoPtr<InjectionModel<CloudType> >
            (
                new MyKinematicLookupTableInjection<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~MyKinematicLookupTableInjection();


    // Member Functions

        //- Set injector locations when mesh is updated
        virtual void updateMesh();

        //- Return the end-of-injection time
        scalar timeEnd() const;

        //- Number of parcels to introduce relative to SOI
        virtual label parcelsToInject(const scalar time0, const scalar time1);

        //- Volume of parcels to introduce relative to SOI
        virtual scalar volumeToInject(const scalar time0, const scalar time1);



        // Injection geometry

            //- Set the injection position and owner cell, tetFace and tetPt
            virtual void setPositionAndCell
            (
                const label parcelI,
                const label nParcels,
                const scalar time,
                vector& position,
                label& cellOwner,
                label& tetFaceI,
                label& tetPtI
            );

            //- Set the parcel properties
            virtual void setProperties
            (
                const label parcelI,
                const label nParcels,
                const scalar time,
                typename CloudType::parcelType& parcel
            );

			//Cholsky decomposition
			tensor CholskyDecompose(symmTensor& sisj);

            //- Flag to identify whether model fully describes the parcel
            virtual bool fullyDescribed() const;
            //{
            //    return true;
            //}

            //- Return flag to identify whether or not injection of parcelI is
            //  permitted
            virtual bool validInjection(const label parcelI);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "MyKinematicLookupTableInjection.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
