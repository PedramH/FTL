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

#include "ParticleStatistics.H"
#include "Pstream.H"
#include "ListListOps.H"
#include "surfaceWriter.H"
#include "Time.H"

// * * * * * * * * * * * * * Protectd Member Functions * * * * * * * * * * * //

template<class CloudType>
void Foam::ParticleStatistics<CloudType>::write()
{
    // check time range
    const fvMesh& mesh = this->owner().mesh();
    if((mesh.time().value() < timeStart_) || mesh.time().value() > timeEnd_)
    {
         return;
    }

    if (ParticleSumPtr_[0].valid()) 
    {
	
		forAll (ParticleSumPtr_, jField)
		{         
		  // write dt*Particle number
		  ParticleSumPtr_[jField]->write();

		  volScalarField& ParticleSum        = ParticleSumPtr_[jField]();
		  scalarField& ndtPart= ParticleSum.internalField(); 

		  volScalarField& ParticleCSum = ParticleCSumPtr_[jField]();

		  ParticleCSum.internalField() = ndtPart / mesh.V();

		  // write concentration field
		  ParticleCSumPtr_[jField]->write();
		  Info << "Writing ParticleStatistics: CMean";
		 
		  if (UMean_)
		  {
	 		  volVectorField& ParticleUSum       = ParticleUSumPtr_[jField]();
			  // create temporary field to store mean U 
			  volVectorField tPartUMean
			  (
				 IOobject
				 (
				     this->modelName() + ":ParticleUMeanG" + name(jField+1),
				     mesh.time().timeName(),
				     mesh,
				     IOobject::NO_READ,
				     IOobject::NO_WRITE
				 ),
				 ParticleUSum
			  );
			  vectorField& PUSum=tPartUMean.internalField();

			  forAll (PUSum, cellI)
			  {
				 scalar ndtp =ndtPart[cellI];
				 if(ndtp > 0) PUSum[cellI]     = PUSum[cellI]/ndtp;
			  }

			  tPartUMean.write();
			  Info << " and UMean" ;

			  if (uPrime2Mean_)
			  {
					volSymmTensorField& ParticleUsqSum     = ParticleUsqSumPtr_[jField]();
					volSymmTensorField tPartUsqMean
					(
						IOobject
						(
						    this->modelName() + ":ParticleUMean2PrimeG" + name(jField+1),
						    mesh.time().timeName(),
						    mesh,
						    IOobject::NO_READ,
						    IOobject::NO_WRITE
						),
						ParticleUsqSum
					);
					symmTensorField& PUsqSum=tPartUsqMean.internalField();

	 				forAll (PUsqSum, cellI)
					{
						 scalar ndtp =ndtPart[cellI];
						 if(ndtp > 0)
						 {
							 vector UMean = PUSum[cellI];
							 PUsqSum[cellI]   = PUsqSum[cellI]/ndtp -
												symmTensor(UMean[0]*UMean[0],UMean[0]*UMean[1],UMean[0]*UMean[2],
														   UMean[1]*UMean[1],UMean[1]*UMean[2],UMean[2]*UMean[2]);
						 }
					}

					tPartUsqMean.write();
					Info << " and UPrime2Mean";

			  }

		  }     

		  //calculate the total number of added particles. Use gSum to make sure to take all the processors 
		  //scalar nPartSum = gSum(ParticleSum);

		  if(writed10_)
		  { 
			  volScalarField& Particled10Sum     = Particled10SumPtr_[jField]();
			  volScalarField tPartd10
			  (
				 IOobject
				 (
				     this->modelName() + ":Particled10G" + name(jField+1),
				     mesh.time().timeName(),
				     mesh,
				     IOobject::NO_READ,
				     IOobject::NO_WRITE
				 ),
				 Particled10Sum
			  );
			  scalarField& Pd10Sum=tPartd10.internalField();

			  forAll (Pd10Sum, cellI)
			  {
				 scalar ndtp =ndtPart[cellI];
				 if(ndtp > 0)
				 {
				   scalar d10     = Pd10Sum[cellI]/ndtp;
				   Pd10Sum[cellI]     = d10;
				 }
			  }

			  tPartd10.write();
			  Info << " and d10";
		  }
		  if(writed32_)
		  { 
			  volScalarField& Particled2Sum      = Particled2SumPtr_[jField]();
		  	  volScalarField& Particled3Sum      = Particled3SumPtr_[jField]();
			  volScalarField tPartd32
			  (
				 IOobject
				 (
				     this->modelName() + ":Particled32G" + name(jField+1),
				     mesh.time().timeName(),
				     mesh,
				     IOobject::NO_READ,
				     IOobject::NO_WRITE
				 ),
				 Particled3Sum
			  );
			  scalarField& Pd32Sum=tPartd32.internalField();

			  scalarField& ndtd2Part= Particled2Sum.internalField();

			  forAll (Pd32Sum, cellI)
			  {
				 scalar ndtd2p =ndtd2Part[cellI];
				 if(ndtd2p > 0)
				 {
				   scalar d32     = Pd32Sum[cellI]/ndtd2p;
				   Pd32Sum[cellI]     = d32;
				 }
			  }

			  tPartd32.write();
			  Particled2SumPtr_[jField]->write();
			  Particled3SumPtr_[jField]->write();
			  Info << " and d32";
		  }
		} //forAll
    }
    else
    {
        FatalErrorIn("void Foam::ParticleStatistics<CloudType>::write()")
            << "ParticleSumPtr not valid" << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ParticleStatistics<CloudType>::ParticleStatistics
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName
)
:
	CloudFunctionObject<CloudType>(dict, owner, modelName, typeName),
/*    ParticleSumPtr_(NULL),
	ParticleCSumPtr_(NULL),
    ParticleUSumPtr_(NULL),
    ParticleUsqSumPtr_(NULL),
    Particled10SumPtr_(NULL),
	Particled2SumPtr_(NULL),
	Particled3SumPtr_(NULL), */
    resetFields_(this->coeffDict().lookup("resetFields")),
	gNum_(readScalar(this->coeffDict().lookup("gNum"))),
	UMean_(this->coeffDict().lookup("UMean")),
    uPrime2Mean_(this->coeffDict().lookup("uPrime2Mean")),
    writed10_(this->coeffDict().lookup("d10")),
	writed32_(this->coeffDict().lookup("d32")),
    timeStart_(readScalar(this->coeffDict().lookup("timeStart"))),
    timeEnd_(readScalar(this->coeffDict().lookup("timeEnd")))
{

	nullifyPtrLists();

	dgMin_.setSize(gNum_);
	dgMax_.setSize(gNum_);
	dgMin_ = readList<scalar>(this->coeffDict().lookup("dgMin"));
	dgMax_ = readList<scalar>(this->coeffDict().lookup("dgMax"));

}


template<class CloudType>
Foam::ParticleStatistics<CloudType>::ParticleStatistics
(
    const ParticleStatistics<CloudType>& ps
)
:
    CloudFunctionObject<CloudType>(ps),
/*    ParticleSumPtr_(NULL),
	ParticleCSumPtr_(NULL),
    ParticleUSumPtr_(NULL),
    ParticleUsqSumPtr_(NULL),
	Particled10SumPtr_(NULL),
	Particled2SumPtr_(NULL),
	Particled3SumPtr_(NULL), */
    resetFields_(ps.resetFields_),
	gNum_(ps.gNum_),
	UMean_(ps.UMean_),
    uPrime2Mean_(ps.uPrime2Mean_),
	writed10_(ps.writed10_),
    writed32_(ps.writed32_),
    timeStart_(ps.timeStart_),
    timeEnd_(ps.timeEnd_),
	dgMin_(ps.dgMin_),
	dgMax_(ps.dgMax_)
{
	nullifyPtrLists();	
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ParticleStatistics<CloudType>::~ParticleStatistics()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::ParticleStatistics<CloudType>::nullifyPtrLists()
{
/*	ParticleSumPtr_.setSize(gNum_,NULL);
	ParticleCSumPtr_.setSize(gNum_,NULL);
	ParticleUSumPtr_.setSize(gNum_,NULL);
	ParticleUsqSumPtr_.setSize(gNum_,NULL);
	Particled10SumPtr_.setSize(gNum_,NULL);
	Particled2SumPtr_.setSize(gNum_,NULL);
	Particled3SumPtr_.setSize(gNum_,NULL);
*/
	ParticleSumPtr_.setSize(gNum_);
	ParticleCSumPtr_.setSize(gNum_);
	ParticleUSumPtr_.setSize(gNum_);
	ParticleUsqSumPtr_.setSize(gNum_);
	Particled10SumPtr_.setSize(gNum_);
	Particled2SumPtr_.setSize(gNum_);
	Particled3SumPtr_.setSize(gNum_);

	forAll (ParticleSumPtr_, jField)
	{
		ParticleSumPtr_[jField].set(NULL);
		ParticleCSumPtr_[jField].set(NULL);
		ParticleUSumPtr_[jField].set(NULL);
		ParticleUsqSumPtr_[jField].set(NULL);
		Particled10SumPtr_[jField].set(NULL);
		Particled2SumPtr_[jField].set(NULL);
		Particled3SumPtr_[jField].set(NULL);
	}
}

template<class CloudType>
void Foam::ParticleStatistics<CloudType>::preEvolve()
{

    if (!ParticleSumPtr_[0].valid())
    {
		const fvMesh& mesh = this->owner().mesh();

		Info << " Paticle statistics is copmuted for Ng = " << gNum_ 
			 <<	"groups: "<< nl;
		Info << " dgMin =" << dgMin_ << nl;
		Info << " dgMax =" << dgMax_ << nl;

		Info << "Creating cloudFunction ParticleStatistics CMean summation array"  << nl;

		forAll (ParticleSumPtr_, jField)
		{
			ParticleSumPtr_[jField].reset
		    (
		        new volScalarField
		        (
		            IOobject
		            (
		                this->modelName() + ":ParticleSumG" + name(jField+1),
		                mesh.time().timeName(mesh.time().startTime().value()),
		                mesh,
		                IOobject::READ_IF_PRESENT,
		                IOobject::NO_WRITE
		            ),
		            mesh,
		            dimensionedScalar("zero", dimless, 0.0)
		        )
		    );
			volScalarField& ParticleSum     = ParticleSumPtr_[jField]();
			scalarField&    nPart=ParticleSum.internalField();
			if(resetFields_)
		    {
		       Info << "Reseting ParticleStatistics averaging fields" << nl;
		       nPart  =scalar(0);
		    }
	
		    ParticleCSumPtr_[jField].reset
		    (
		        new volScalarField
		        (
		            IOobject
		            (
		                this->modelName() + ":ParticleCMeanG" + name(jField+1),
		                mesh.time().timeName(mesh.time().startTime().value()),
		                mesh,
		                IOobject::NO_READ,
		                IOobject::NO_WRITE
		            ),
		            mesh,
		            dimensionedScalar("zero", dimless/dimVolume, 0.0)
		        )
		    );
			volScalarField& ParticleCSum    = ParticleCSumPtr_[jField]();
			scalarField&    PCSum=ParticleCSum.internalField();

			if(resetFields_) PCSum  =scalar(0);

			if(UMean_)
			{
				Info << "Creating cloudFunction ParticleStatistics UMean summation array"  << nl;
				// since the cloudfunction is called in the first iteration, mesh.time() referces
				// not to the start time, but to the first next time step. You want to have the start time value in case
				// you want to continue with the previous run, therefore I enfore to take mesh.time.startTime()
				// This is not string (as timeName) but a dimensionedscalar. Therefor you take the value and pass
				// it to the Time.H function timeName to turn the scalar into a string with
				// mesh.time().timeName(mesh.time().startTime().value())
				//
				ParticleUSumPtr_[jField].reset
				(
					new volVectorField
					(
					    IOobject
					    (
					        this->modelName() + ":ParticleUMeanG" + name(jField+1),
					        mesh.time().timeName(mesh.time().startTime().value()),
					        mesh,
					        IOobject::READ_IF_PRESENT,
					        IOobject::NO_WRITE
					    ),
					    mesh,
					    dimensionedVector("zero", dimVelocity, vector(0.0,0,0))
					)
				);
				volVectorField& ParticleUSum    = ParticleUSumPtr_[jField]();
				vectorField&    PUSum=ParticleUSum.internalField();

				if(uPrime2Mean_)
				{
					Info << "Creating cloudFunction ParticleStatistics UPrime2Mean summation array"  << nl;
					ParticleUsqSumPtr_[jField].reset
					(
						new volSymmTensorField
						(
							IOobject
							(
							    this->modelName() + ":ParticleUMean2PrimeG" + name(jField+1),
							    mesh.time().timeName(mesh.time().startTime().value()),
							    mesh,
							    IOobject::READ_IF_PRESENT,
							    IOobject::NO_WRITE
							),
							mesh,
							dimensionedSymmTensor("zero", dimVelocity*dimVelocity, symmTensor(0,0,0,0,0,0))
						)
					);
					volSymmTensorField& ParticleUsqSum  = ParticleUsqSumPtr_[jField]();
					symmTensorField&    PUsqSum=ParticleUsqSum.internalField();

					if(resetFields_)
					{
					   PUSum  =vector(0,0,0);
					   PUsqSum=symmTensor(0,0,0,0,0,0);
					}
					else
					{

					  // the average fields are read, but the contain the averages <U>=Sum U_i/Sum n_i  and 
					  // <u'2> = Sum_i U_i*U_i/Sum n_i - <U>^2 (where U=<U>+u')
					  forAll (PUSum, cellI)
					  {
						scalar np       = nPart[cellI];
						vector UMean    = PUSum[cellI];
						PUsqSum[cellI]  = np*(PUsqSum[cellI]+symmTensor(UMean[0]*UMean[0],UMean[0]*UMean[1],
								UMean[0]*UMean[2],UMean[1]*UMean[1],UMean[1]*UMean[2],UMean[2]*UMean[2]));
						PUSum[cellI]       = np*UMean;
					  }
					}
				}
				else
				{
					if(resetFields_)
					{
					   PUSum  =vector(0,0,0);
					}
					else
					{
					  // the average fields are read, but the contain the averages <U>=Sum U_i/Sum n_i  and 
					  // <u'2> = Sum_i U_i*U_i/Sum n_i - <U>^2 (where U=<U>+u')
					  forAll (PUSum, cellI)
					  {
						scalar np       = nPart[cellI];
						PUSum[cellI]       = np*PUSum[cellI];
					  }
					}
				}

			} //if(UMean_)

			if(writed10_)
			{
				Info << "Creating cloudFunction ParticleStatistics d10 summation array"  << nl;
				Particled10SumPtr_[jField].reset
				(
				    new volScalarField
				    (
				        IOobject
				        (
				            this->modelName() + ":Particled10G" + name(jField+1),
				            mesh.time().timeName(mesh.time().startTime().value()),
				            mesh,
				            IOobject::READ_IF_PRESENT,
				            IOobject::NO_WRITE
				        ),
				        mesh,
				        dimensionedScalar("zero", dimLength, 0.0)
				    )
				);
				volScalarField& Particled10Sum    = Particled10SumPtr_[jField]();
				scalarField&    Pd10Sum=Particled10Sum.internalField();

				if(resetFields_)
				{
				   Pd10Sum=scalar(0);
				}
				else
				{
				  forAll (Pd10Sum, cellI)
				  {
				    scalar np       = nPart[cellI];
					Pd10Sum[cellI]       = np*Pd10Sum[cellI];
				  }
				}
			}
		
			if(writed32_)
			{
				Info << "Creating cloudFunction ParticleStatistics d32 summation array"  << nl;
				Particled2SumPtr_[jField].reset
				(
				    new volScalarField
				    (
				        IOobject
				        (
				            this->modelName() + ":Particled2G" + name(jField+1),
				            mesh.time().timeName(mesh.time().startTime().value()),
				            mesh,
				            IOobject::READ_IF_PRESENT,
				            IOobject::NO_WRITE
				        ),
				        mesh,
				        dimensionedScalar("zero", dimArea, 0.0)
				    )
				);
				volScalarField& Particled2Sum    = Particled2SumPtr_[jField]();
				scalarField&    Pd2Sum=Particled2Sum.internalField();

				Particled3SumPtr_[jField].reset
				(
				    new volScalarField
				    (
				        IOobject
				        (
				            this->modelName() + ":Particled3G" + name(jField+1),
				            mesh.time().timeName(mesh.time().startTime().value()),
				            mesh,
				            IOobject::READ_IF_PRESENT,
				            IOobject::NO_WRITE
				        ),
				        mesh,
				        dimensionedScalar("zero", dimVolume, 0.0)
				    )
				);
				volScalarField& Particled3Sum    = Particled3SumPtr_[jField]();  
				scalarField&    Pd3Sum=Particled3Sum.internalField();

				if(resetFields_)
				{
				   Pd2Sum=scalar(0);
				   Pd3Sum=scalar(0);
				}
			}

		} //forAll
   }

}


template<class CloudType>
void Foam::ParticleStatistics<CloudType>::postEvolve()
{

    CloudFunctionObject<CloudType>::postEvolve();

}


template<class CloudType>
void Foam::ParticleStatistics<CloudType>::postMove
(
    const parcelType& p,
    const label cellI,
    const scalar dt,
	const point& position0,
    bool&
)
{
    const fvMesh& mesh = this->owner().mesh();

    // question: is there at better way to filter on the time range: not in postMove but before
    // you call this function. Now you will always enter this function even if it is not necessary
    if((mesh.time().value()<timeStart_) || mesh.time().value() > timeEnd_)
    {
      // not in the time range: go back
      return;
    }

	label jField;

	for (jField = 0; jField < gNum_; jField++)
	{ 
		if (p.d() > dgMin_[jField] && p.d() <= dgMax_[jField]) 
		{
			break;
		}
	}

	if (jField == gNum_) 
	{
      // not in the specified gropus: go back
      return;
    }

	scalar partN = dt/mesh.time().deltaTValue()*p.nParticle();

	volScalarField& ParticleSum     = ParticleSumPtr_[jField]();

    ParticleSum[cellI]  += partN; //p.nParticle();
	if(UMean_) 
	{
		volVectorField& ParticleUSum    = ParticleUSumPtr_[jField]();
		ParticleUSum[cellI] += partN*p.U();
		if(uPrime2Mean_)
		{
		   volSymmTensorField& ParticleUsqSum  = ParticleUsqSumPtr_[jField]();
		  // keep the sum of the Ux^2 , Uy^2 and Uz^2, ... in a vector
		  ParticleUsqSum[cellI] += partN*
		    symmTensor(
		           p.U()[0]*p.U()[0],
				   p.U()[0]*p.U()[1],
				   p.U()[0]*p.U()[2],
		           p.U()[1]*p.U()[1],
				   p.U()[1]*p.U()[2],
		           p.U()[2]*p.U()[2]
		          );
		}
	}
	if(writed10_)
	{
		volScalarField& Particled10Sum    = Particled10SumPtr_[jField]();
		Particled10Sum[cellI] += partN*p.d();
	}
	if(writed32_) 
	{
		volScalarField& Particled2Sum    = Particled2SumPtr_[jField]();
		volScalarField& Particled3Sum    = Particled3SumPtr_[jField]();
		scalar nd2 = partN*p.d()*p.d();
		Particled2Sum[cellI] += nd2;
		Particled3Sum[cellI] += nd2*p.d();
	}
}


// ************************************************************************* //
