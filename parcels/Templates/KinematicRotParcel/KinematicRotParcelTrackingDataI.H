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

\*---------------------------------------------------------------------------*/

template<class ParcelType>
template<class CloudType>
inline Foam::KinematicRotParcel<ParcelType>::TrackingData<CloudType>::TrackingData
(
    CloudType& cloud,
    trackPart part,
	fieldUpdate updates //added1
)
:
    ParcelType::template TrackingData<CloudType>(cloud, part, updates), //added1
	//particleStatistics2_(false), //added1
	OmegaPStat_(NULL) 
	/*//OmegaC memory
    curlUInterp_
    (
        interpolation<vector>::New
        (
            cloud.solution().interpolationSchemes(),
            cloud.curlU()
        )
    )*/
{
	//if (cloud.stat2Comp()) particleStatistics2_ = true; //added1
	
	if (cloud.preMotStat2Comp() || cloud.preColStat2Comp())
	{
		OmegaPStat_.reset
		(
			new DimensionedField<vector, volMesh>
			(
			    IOobject
			    (
			        cloud.name() + ":OmegaPStat",
			        cloud.db().time().timeName(),
			        cloud.db(),
			        IOobject::NO_READ,
	            	IOobject::NO_WRITE
			    ),
			    cloud.mesh(),
			    dimensionedVector("zero", dimless/dimTime, vector::zero)
			)
		);
	}

	/*if (particleStatistics1_)
    {	
		uAverage_.reset
		(
		    AveragingMethod<vector>::New
		    (
		        IOobject
		        (
		            cloud.name() + ":uAverage",
		            cloud.db().time().timeName(),
		            cloud.mesh()
		        ),
		        cloud.solution().dict(),
		        cloud.mesh()
		    )
		);

		if (particleStatistics2_)
    	{
			u2Average_.reset
			(
				AveragingMethod<vector>::New
				(
				    IOobject
				    (
				        cloud.name() + ":u2Average",
				        cloud.db().time().timeName(),
				        cloud.mesh()
				    ),
				    cloud.solution().dict(),
				    cloud.mesh()
				)
			);
			dAverage_.reset
			(
				AveragingMethod<scalar>::New
				(
				    IOobject
				    (
				        cloud.name() + ":dAverage",
				        cloud.db().time().timeName(),
				        cloud.mesh()
				    ),
				    cloud.solution().dict(),
				    cloud.mesh()
				)
			);
			d2Average_.reset
			(
				AveragingMethod<scalar>::New
				(
				    IOobject
				    (
				        cloud.name() + ":d2Average",
				        cloud.db().time().timeName(),
				        cloud.mesh()
				    ),
				    cloud.solution().dict(),
				    cloud.mesh()
				)
			);
			nAverage_.reset
			(
				AveragingMethod<scalar>::New
				(
				    IOobject
				    (
				        cloud.name() + ":nAverage",
				        cloud.db().time().timeName(),
				        cloud.mesh()
				    ),
				    cloud.solution().dict(),
				    cloud.mesh()
				)
			);
		}
	}*/
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ParcelType>
template<class CloudType>
Foam::KinematicRotParcel<ParcelType>::TrackingData<CloudType>::~TrackingData()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//OmegaC memory
/*template<class ParcelType>
template<class CloudType>
inline const Foam::interpolation<Foam::vector>&
Foam::KinematicRotParcelParcel<ParcelType>::TrackingData<CloudType>::curlUInterp() const
{
    return curlUInterp_();
}*/


template<class ParcelType> //added1
template<class CloudType>
inline void
Foam::KinematicRotParcel<ParcelType>::TrackingData<CloudType>::updateParticleStatistics
(
    //CloudType& cloud
)
{
	ParcelType::template TrackingData<CloudType>::updateParticleStatistics(); //(cloud);
	
	CloudType& cloud = this->cloud();
	if (this->particleStatistics2_)
    {    
		Info<< "For step:" << this->updates() << "; OmegaPStat is updated." << nl; //updates() does not work
		
		vectorField& OmegaP = OmegaPStat_().field();

		OmegaP = vector::zero;

		autoPtr<DimensionedField<scalar, volMesh> > IPPtr
		(
			new DimensionedField<scalar, volMesh>
			(
			    IOobject
			    (
			        cloud.name() + ":IP",
			        cloud.db().time().timeName(),
			        cloud.db(),
			        IOobject::NO_READ,
	            	IOobject::NO_WRITE
			    ),
			    cloud.mesh(),
			    dimensionedScalar("zero", dimMass*dimLength*dimLength/dimVolume, 0.0)
			)
		);
		//DimensionedField<scalar, volMesh>& MP = MPPtr();
		scalarField& IP = IPPtr().field();

		// averaging sums Sommerfeld's collision model
		IP = 0;
		forAllConstIter(typename CloudType, cloud, iter)
		{
			const typename CloudType::parcelType& p = iter();
			const label cellI = p.cell();
			//const tetIndices tetIs(cellI, p.tetFace(), p.tetPt(), cloud.mesh());

			const scalar n = p.nParticle();
			const scalar I = n*p.momentOfInertia();

			IP[cellI] += I;
			OmegaP[cellI] += I*p.Omega();
		}
		OmegaP /= max(IP, SMALL);
	}
}


/*
template<class ParcelType>
template<class CloudType>
inline void
Foam::KinematicRotParcel<ParcelType>::TrackingData<CloudType>::updateAverages
(
    CloudType& cloud
)
{
	if (particleStatistics1_)
    {    
		// zero the sums
		uAverage_() = vector::zero;

		// temporary weights
		autoPtr<AveragingMethod<scalar> > weightAveragePtr
		(
		    AveragingMethod<scalar>::New
		    (
		        IOobject
		        (
		            cloud.name() + ":weightAverage",
		            cloud.db().time().timeName(),
		            cloud.mesh()
		        ),
		        cloud.solution().dict(),
		        cloud.mesh()
		    )
		);
		AveragingMethod<scalar>& weightAverage = weightAveragePtr();

		// averaging sums //CRW dispersion and Sommerfeld's collision model
		weightAverage = 0;
		forAllConstIter(typename CloudType, cloud, iter)
		{
		    const typename CloudType::parcelType& p = iter();
		    const tetIndices tetIs(p.cell(), p.tetFace(), p.tetPt(), cloud.mesh());

		    const scalar m = p.nParticle()*p.mass();

		    weightAverage.add(p.position(), tetIs, m);
		    uAverage_->add(p.position(), tetIs, m*p.U());
		}
		weightAverage.average();
		uAverage_->average(weightAverage);

		if (particleStatistics2_)
   		{ 

			u2Average_() = vector::zero;
			dAverage_() = 0;
			d2Average_() = 0;
			nAverage_() = 0;


			// Sommerfeld's collision model
			forAllConstIter(typename CloudType, cloud, iter)
			{
				const typename CloudType::parcelType& p = iter();
				const tetIndices tetIs(p.cell(), p.tetFace(), p.tetPt(), cloud.mesh());

				const scalar m = p.nParticle()*p.mass();

				weightAverage.add(p.position(), tetIs, m);
				u2Average_->add(p.position(), tetIs, m*p.U()*p.U());
		
				nAverage_->add(p.position(), tetIs, p.nParticle());
				dAverage_->add(p.position(), tetIs, p.nParticle()*p.d());
				d2Average_->add(p.position(), tetIs, p.nParticle()*p.d()*p.d());		
			}
			u2Average_->average(weightAverage);
			nAverage_->average();
			dAverage_->average(nAverage_);
			d2Average_->average(nAverage_);
		}
	}
}
*/


// ************************************************************************* //
