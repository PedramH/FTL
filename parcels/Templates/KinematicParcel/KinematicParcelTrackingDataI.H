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

\*---------------------------------------------------------------------------*/

template<class ParcelType>
template<class CloudType>
inline Foam::KinematicParcel<ParcelType>::TrackingData<CloudType>::TrackingData
(
    CloudType& cloud,
    trackPart part,
	fieldUpdate updates //added1
)
:
    ParcelType::template TrackingData<CloudType>(cloud),
    rhoInterp_
    (
        interpolation<scalar>::New
        (
            cloud.solution().interpolationSchemes(),
            cloud.rho()
        )
    ),
    UInterp_
    (
        interpolation<vector>::New
        (
            cloud.solution().interpolationSchemes(),
            cloud.U()
        )
    ),
    muInterp_
    (
        interpolation<scalar>::New
        (
            cloud.solution().interpolationSchemes(),
            cloud.mu()
        )
    ),
    g_(cloud.g().value()),
    part_(part),
	updates_(updates), //added1
	particleStatistics1_(false), 
	particleStatistics2_(false), //added1
	UPStat_(NULL),
	uPuPStat_(NULL), //added1
	NPStat_(NULL),
	DPStat_(NULL),
	dPdPStat_(NULL)
{
	//if (cloud.stat1Comp()) particleStatistics1_ = true; //added1

	if (cloud.preMotStat1Comp() || cloud.preColStat1Comp()) //added1
    {
		UPStat_.reset
		(
		    new DimensionedField<vector, volMesh>
		    (
		        IOobject
		        (
		            cloud.name() + ":UPStat",
		            cloud.db().time().timeName(),
		            cloud.db(),
		            IOobject::NO_READ,
                	IOobject::NO_WRITE
		        ),
		        cloud.mesh(),
		        dimensionedVector("zero", dimVelocity, vector::zero)
		    )
		);
	}

	if (cloud.preMotStat2Comp() || cloud.preColStat2Comp())
	{
		uPuPStat_.reset
		(
			new DimensionedField<vector, volMesh>
			(
			    IOobject
			    (
			        cloud.name() + ":uPuPStat",
			        cloud.db().time().timeName(),
			        cloud.db(),
			        IOobject::NO_READ,
	            	IOobject::NO_WRITE
			    ),
			    cloud.mesh(),
			    dimensionedVector("zero", dimVelocity*dimVelocity, vector::zero)
			)
		);

		NPStat_.reset
		(
			new DimensionedField<scalar, volMesh>
			(
			    IOobject
			    (
			        cloud.name() + ":NPStat",
			        cloud.db().time().timeName(),
			        cloud.db(),
			        IOobject::NO_READ,
	            	IOobject::NO_WRITE
			    ),
			    cloud.mesh(),
			    dimensionedScalar("zero", dimless/dimVolume, 0.0)
			)
		);

		DPStat_.reset
		(
			new DimensionedField<scalar, volMesh>
			(
			    IOobject
			    (
			        cloud.name() + ":DPStat",
			        cloud.db().time().timeName(),
			        cloud.db(),
			        IOobject::NO_READ,
	            	IOobject::NO_WRITE
			    ),
			    cloud.mesh(),
			    dimensionedScalar("zero", dimLength, 0.0)
			)
		);

		dPdPStat_.reset
		(
			new DimensionedField<scalar, volMesh>
			(
			    IOobject
			    (
			        cloud.name() + ":dPdPStat",
			        cloud.db().time().timeName(),
			        cloud.db(),
			        IOobject::NO_READ,
	            	IOobject::NO_WRITE
			    ),
			    cloud.mesh(),
			    dimensionedScalar("zero", dimLength*dimLength, 0.0)
			)
		);
	}
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ParcelType>
template<class CloudType>
Foam::KinematicParcel<ParcelType>::TrackingData<CloudType>::~TrackingData()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParcelType>
template<class CloudType>
inline const Foam::interpolation<Foam::scalar>&
Foam::KinematicParcel<ParcelType>::TrackingData<CloudType>::rhoInterp() const
{
    return rhoInterp_();
}


template<class ParcelType>
template<class CloudType>
inline const Foam::interpolation<Foam::vector>&
Foam::KinematicParcel<ParcelType>::TrackingData<CloudType>::UInterp() const
{
    return UInterp_();
}


template<class ParcelType>
template<class CloudType>
inline const Foam::interpolation<Foam::scalar>&
Foam::KinematicParcel<ParcelType>::TrackingData<CloudType>::muInterp() const
{
    return muInterp_();
}


template<class ParcelType>
template<class CloudType>
inline const Foam::vector&
Foam::KinematicParcel<ParcelType>::TrackingData<CloudType>::g() const
{
    return g_;
}


template<class ParcelType>
template<class CloudType>
inline typename Foam::KinematicParcel<ParcelType>::template
TrackingData<CloudType>::trackPart
Foam::KinematicParcel<ParcelType>::TrackingData<CloudType>::part() const
{
    return part_;
}


template<class ParcelType>
template<class CloudType>
inline typename Foam::KinematicParcel<ParcelType>::template
TrackingData<CloudType>::trackPart&
Foam::KinematicParcel<ParcelType>::TrackingData<CloudType>::part()
{
    return part_;
}


template<class ParcelType> //added1
template<class CloudType>
inline typename Foam::KinematicParcel<ParcelType>::template
TrackingData<CloudType>::fieldUpdate
Foam::KinematicParcel<ParcelType>::TrackingData<CloudType>::updates() const
{
    return updates_;
}


template<class ParcelType> //added1
template<class CloudType>
inline typename Foam::KinematicParcel<ParcelType>::template
TrackingData<CloudType>::fieldUpdate&
Foam::KinematicParcel<ParcelType>::TrackingData<CloudType>::updates()
{
    return updates_;
}


template<class ParcelType> //added1
template<class CloudType>
inline void
Foam::KinematicParcel<ParcelType>::TrackingData<CloudType>::updateParticleStatistics
(
    //CloudType& cloud
)
{
	CloudType& cloud = this->cloud();
	particleStatistics1_ = false;
	particleStatistics2_ = false;
	switch (updates_)
    {
        case motion:
        {
            particleStatistics1_ = cloud.preMotStat1Comp();
			particleStatistics2_ = cloud.preMotStat2Comp();

            break;
        }

        case statCollision:
        {
            particleStatistics1_ = cloud.preColStat1Comp();
			particleStatistics2_ = cloud.preColStat2Comp();

            break;
        }

        default:
        {
            FatalErrorIn
            (
                "KinematicParcel<ParcelType>::updateParticleStatistics(CloudType& cloud)"
            )   << updates_ << " is an invalid value for field updates."
                << abort(FatalError);
        }
    }


	if (particleStatistics1_)
    {    
		Info<< "For step:" << updates_ << "; UPStat is updated." << nl;
		// zero the sums
		vectorField& UP = UPStat_().field();
		UP = vector::zero;

		// temporary weights
		autoPtr<DimensionedField<scalar, volMesh> > MPPtr
		(
			new DimensionedField<scalar, volMesh>
			(
			    IOobject
			    (
			        cloud.name() + ":MP",
			        cloud.db().time().timeName(),
			        cloud.db(),
			        IOobject::NO_READ,
	            	IOobject::NO_WRITE
			    ),
			    cloud.mesh(),
			    dimensionedScalar("zero", dimMass, 0.0)
			)
		);
		//DimensionedField<scalar, volMesh>& MP = MPPtr();
		scalarField& MP = MPPtr().field();

		// averaging sums //CRW dispersion and Sommerfeld's collision model
		MP = 0;
		forAllConstIter(typename CloudType, cloud, iter)
		{
		    const typename CloudType::parcelType& p = iter();
			const label cellI = p.cell();
		    //const tetIndices tetIs(cellI, p.tetFace(), p.tetPt(), cloud.mesh());

		    const scalar m = p.nParticle()*p.mass();

			MP[cellI] += m;
		    UP[cellI] += m*p.U();
		}
		UP /= max(MP, SMALL);

		if (particleStatistics2_)
		{    
			Info<< "For step:" << updates_ << "; uPuPStat, NPStat, DPStat and dPdPStat are updated." << nl;
		
			vectorField& uPuP = uPuPStat_().field();
			scalarField& NP = NPStat_().field();
			scalarField& DP = DPStat_().field();
			scalarField& dPdP = dPdPStat_().field();

			uPuP = vector::zero;
			NP = 0;
			DP = 0;
			dPdP = 0;

			// averaging sums //CRW dispersion and Sommerfeld's collision model
			forAllConstIter(typename CloudType, cloud, iter)
			{
				const typename CloudType::parcelType& p = iter();
				const label cellI = p.cell();
				//const tetIndices tetIs(cellI, p.tetFace(), p.tetPt(), cloud.mesh());

				const scalar n = p.nParticle();
				const scalar m = n*p.mass();
			
				uPuP[cellI] += m*cmptMultiply(p.U(),p.U());	
				NP[cellI] += n;
				DP[cellI] += n*p.d();
				dPdP[cellI] += n*p.d()*p.d();
			}
			uPuP /= max(MP, SMALL); //UP2
			uPuP = uPuP - cmptMultiply(UP,UP); //standard deviation
			DP /= max(NP, SMALL); 
			dPdP /= max(NP, SMALL); 
			dPdP = dPdP - DP*DP; //standard deviation
			NP /= cloud.mesh().V();
		}
	}
}


// ************************************************************************* //
