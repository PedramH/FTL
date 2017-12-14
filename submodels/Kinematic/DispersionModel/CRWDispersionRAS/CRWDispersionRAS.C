/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "CRWDispersionRAS.H"
#include "demandDrivenData.H"
#include "turbulenceModel.H" //changed
#include "fvc.H"
//#include "Scalar.H"
#include "myAxesRotation.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class CloudType>
Foam::tmp<Foam::volSymmTensorField>
Foam::CRWDispersionRAS<CloudType>::RfModel()
{
	//changed
    const objectRegistry& obr = this->owner().mesh();
    const word turbName =
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            this->owner().U().group()
        );

    if (obr.foundObject<turbulenceModel>(turbName))
    {
        const turbulenceModel& model =
            obr.lookupObject<turbulenceModel>(turbName);
        return model.R();
    }
    else
    {
        FatalErrorIn
        (
            "Foam::tmp<Foam::volScalarField>"
            "Foam::DispersionRASModel<CloudType>::RfModel() const"
        )
            << "Turbulence model not found in mesh database" << nl
            << "Database objects include: " << obr.sortedToc()
            << abort(FatalError);

        return tmp<volSymmTensorField>(NULL);
    }
}


template<class CloudType>
Foam::tmp<Foam::volVectorField>
Foam::CRWDispersionRAS<CloudType>::driftModel()
{
//Info << "M2 ownRf_" << ownRf_ << nl;
		if (ownRf_ && anisotropy_) //????
        {
			return tmp<volVectorField>
					(
						new volVectorField
						(
							IOobject
							(
								"dispersionDrift",
							    RfPtr_->mesh().time().timeName(),
							    RfPtr_->mesh(),
							    IOobject::NO_READ,
				 				IOobject::NO_WRITE
							),
							fvc::div(*RfPtr_) 
						)
					);
        }
        else
        {
            return tmp<volVectorField>
					(
						new volVectorField
						(
							IOobject
							(
								"dispersionDrift",
							    this->kPtr_->mesh().time().timeName(),
							    this->kPtr_->mesh(),
							    IOobject::NO_READ,
				 				IOobject::NO_WRITE
							),
							2.0/3*fvc::grad(*this->kPtr_)
						)
					);
        }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::CRWDispersionRAS<CloudType>::CRWDispersionRAS
(
    const dictionary& dict,
    CloudType& owner
)
:
	DispersionRASModel<CloudType>(dict, owner, typeName), //DispersionRASModel<CloudType>(dict, owner),
	ranGen_(label(0)),
	CL_(this->coeffDict().lookupOrDefault("CL",0.4819)),
	beta_(this->coeffDict().lookupOrDefault("beta",0.356)),
	anisotropy_(this->coeffDict().lookupOrDefault("anisotropy",false)),
	uiuiMin_(this->coeffDict().lookupOrDefault("uiuiMin",ROOTVSMALL)),
	nearWallT_(this->coeffDict().lookupOrDefault("nearWallT",false)),
	Cmu_(this->coeffDict().lookupOrDefault("Cmu",0.09)),
	yPlusBL_(this->coeffDict().lookupOrDefault("yPlusBL",30.0)),
	RfPtr_(NULL),
    ownRf_(false),
	driftPtr_(NULL),
    ownDrift_(false),
	UpMean_(NULL),
	yr_(this->owner().mesh()),
	kInterp_(NULL),
	epsInterp_(NULL),
	RfInterp_(NULL),
	driftInterp_(NULL),
	nwInterp_
    (
        interpolation<vector>::New
        (
            this->owner().solution().interpolationSchemes(),
            yr_.n()
        )
    ),
	yInterp_
    (
        interpolation<scalar>::New
        (
            this->owner().solution().interpolationSchemes(),
            yr_.y()
        )
    )
{
	CL_=1.5*CL_; //aniso based code

	//change cloud flag to compute statistics group 1
	owner.preMotStat1Comp() = true;
}


template<class CloudType>
Foam::CRWDispersionRAS<CloudType>::CRWDispersionRAS
(
    CRWDispersionRAS<CloudType>& dm
)
:
    DispersionRASModel<CloudType>(dm),
	ranGen_(dm.ranGen_),
	CL_(dm.CL_),
	beta_(dm.beta_),
	anisotropy_(dm.anisotropy_),
	uiuiMin_(dm.uiuiMin_),
	nearWallT_(dm.nearWallT_),
	Cmu_(dm.Cmu_),
	yPlusBL_(dm.yPlusBL_),
	RfPtr_(dm.RfPtr_),
    ownRf_(dm.ownRf_),
	driftPtr_(dm.driftPtr_),
    ownDrift_(dm.ownDrift_),
	UpMean_(dm.UpMean_),
	yr_(dm.owner().mesh()),
	nwInterp_(dm.nwInterp_),
	yInterp_(dm.yInterp_)
{
	dm.ownRf_ = false;
	dm.ownDrift_ = false;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::CRWDispersionRAS<CloudType>::~CRWDispersionRAS()
{
	cacheFields(false);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::CRWDispersionRAS<CloudType>::cacheFields(const bool store)
{
	DispersionRASModel<CloudType>::cacheFields(store);

	if (store)
    {
		tmp<volSymmTensorField> tRf = RfModel();
		if (tRf.isTmp())
		{
		    RfPtr_ = tRf.ptr();
		    ownRf_ = true;
		}
		else
		{
		    RfPtr_ = tRf.operator->();
		    ownRf_ = false;
		}
		driftPtr_ = driftModel().ptr();
        ownDrift_ = true;

		//kInterp_->psi() = *this->kPtr_;
		kInterp_ = interpolation<scalar>::New
					(
						this->owner().solution().interpolationSchemes(),
						*this->kPtr_ 
					);
		//epsInterp_->psi() = *this->epsilonPtr_;
		epsInterp_ = interpolation<scalar>::New
					(
						this->owner().solution().interpolationSchemes(),
						*this->epsilonPtr_ 
					);

		//RfInterp_->psi() = *this->RfPtr_;
		RfInterp_ = interpolation<symmTensor>::New
						(
							this->owner().solution().interpolationSchemes(),
							*this->RfPtr_
						);	

		//driftInterp_->psi() = *this->driftPtr_;
		driftInterp_ = interpolation<vector>::New
						(
							this->owner().solution().interpolationSchemes(),
							*this->driftPtr_
						);		

		const fvMesh& mesh = this->owner().mesh();
        const word& cloudName = this->owner().name();

        const DimensionedField<vector, volMesh>& UpMean =
            mesh.lookupObject<DimensionedField<vector, volMesh> >
            (
                cloudName + ":UPStat"
            );
		UpMean_ = &UpMean;

		if (nearWallT_ && this->owner().mesh().changing())
		{
		    yr_.correct();
		} 
	}
	else
	{
		if (ownRf_ && RfPtr_)
        {
            deleteDemandDrivenData(RfPtr_);
			RfPtr_ = NULL;
            ownRf_ = false;
        }
		if (ownDrift_)
        {
            deleteDemandDrivenData(driftPtr_);
            driftPtr_ = NULL;
            ownDrift_ = false;
        }
		UpMean_ = NULL;
	}
}

/*
template<class CloudType>
Foam::vector Foam::CRWDispersionRAS<CloudType>::update
(
    const scalar dt,
    const label cellI,
    const vector& U,
    const vector& Uc,
    vector& UTurb,
    scalar& tTurb
)
{
    const scalar k = this->kPtr_->internalField()[cellI];
    const scalar epsilon =
        this->epsilonPtr_->internalField()[cellI] + ROOTVSMALL;
	vector driftCorr  = this->driftPtr_->internalField()[cellI];

	//Info << "M1 k =" << k << "eps =" << epsilon << "driftCorr =" << driftCorr << nl;

	tTurb = CL_/1.5 * (k+ ROOTVSMALL) / epsilon;

	CRW_Dispersion
	(
		dt,
		cellI,
		U,
		Uc,
		UTurb,
		tTurb,
		k,
		epsilon,
		driftCorr
    );

    return Uc + UTurb;
}
*/
template<class CloudType>
Foam::vector Foam::CRWDispersionRAS<CloudType>::update
(
    const scalar dt,
    const label cellI,
    const vector& U,
    const vector& Uc,
    vector& UTurb,
    scalar& tTurb, //=taop here
	const scalar& muc,
	const scalar& rhoc,
	const typename CloudType::parcelType& p
)
{
	tetIndices tetIs = p.currentTetIndices();

    scalar k = max( kInterp_->interpolate(p.position(), tetIs), ROOTVSMALL);

    scalar epsilon = max( epsInterp_->interpolate(p.position(), tetIs), ROOTVSMALL);

	symmTensor uiujf(symmTensor::zero);

	if (ownRf_ && anisotropy_) //????
	{
		//anisotropic
		uiujf = RfInterp_->interpolate(p.position(), tetIs);
	}
	else
	{
		//isotropic
		uiujf = symmTensor(2*k/3,0,0,2*k/3,0,2*k/3);
	}

	if (nearWallT_)
	{
		const scalar y = max( yInterp_->interpolate(p.position(), tetIs), 0);
		scalar uStar = pow(Cmu_, 0.25) * sqrt(k);
		const scalar yPlus = uStar * y / (muc/rhoc);

		if (yPlus < yPlusBL_) 
		{
			uiujf[0] = sqr(0.4*yPlus/(1+0.0239*pow(yPlus, 1.496)) *uStar);
			uiujf[1] = 0.0;
			uiujf[2] = 0.0;
			uiujf[3] = sqr(0.0116*yPlus*yPlus/(1+0.203*yPlus+0.0014*pow(yPlus, 2.421)) *uStar);
			uiujf[4] = 0.0;
			uiujf[5] = sqr(0.19*yPlus/(1+0.0361*pow(yPlus, 1.322)) *uStar);

			epsilon = 1.0/(4.529+0.0116*pow(yPlus, 1.75)+0.768*sqrt(yPlus));

			const vector nw = nwInterp_->interpolate(p.position(), tetIs);
			vector nt = Uc - (nw & Uc)*nw;
			//nt /= mag(nt) + ROOTVSMALL;

			if (mag(nt)*mag(nw) != 0)
			{
				myAxesRotation AR(nt, nw, 0);
				uiujf = AR.transformTensor(uiujf); //from local to global
			}
		}
	}

	uiujf[0] = max(uiujf[0], ROOTVSMALL);
	uiujf[3] = max(uiujf[3], ROOTVSMALL);
	uiujf[5] = max(uiujf[5], ROOTVSMALL);

	vector driftCorr  = driftInterp_->interpolate(p.position(), tetIs);

	vector UpMean = UpMean_->field()[cellI];

	CRW_Dispersion
	(
		dt,
		cellI,
		U,
		Uc,
		UTurb,
		tTurb, //=taop here
		k,
		epsilon,
		uiujf,
		driftCorr,
		UpMean
    );

    return Uc + UTurb;
}

template<class CloudType>
void Foam::CRWDispersionRAS<CloudType>::tTurbUpdate
(
	scalar& tTurb,
	scalar& taop
)
{
	tTurb = taop;
}

template<class CloudType>
void Foam::CRWDispersionRAS<CloudType>::CRW_Dispersion
(
	const scalar dt,
	const label cellI,
	const vector& U,
	const vector& Uc,
	vector& UTurb,
	scalar& taop,
	const scalar& k,
	const scalar& epsilon,
	symmTensor uiujf,
	vector driftCorr,
	const vector UpMean
)
{
	//cachedRandom& rnd = this->owner().rndGen();

	//vector Up = U;
	vector Urel = UpMean - Uc; //vector Urel = Up - Uc;
	scalar UrelMag = mag(Urel);

			vector randv;
			randv.x() = ranGen_.GaussNormal();
			randv.y() = ranGen_.GaussNormal();
			randv.z() = ranGen_.GaussNormal();

			/*for (int i=0; i<3; i++)
			{
				// Numerical Recipes... Ch. 7. Random Numbers...
		        scalar x1 = 0.0;
		        scalar x2 = 0.0;
		        scalar rsq = 10.0;
		        while ((rsq > 1.0) || (rsq == 0.0))
		        {
		            x1 = 2.0*rnd.sample01<scalar>() - 1.0;
		            x2 = 2.0*rnd.sample01<scalar>() - 1.0;
		            rsq = x1*x1 + x2*x2;
		        }
				randv[i] = rsq;
			}*/

			//rotation of coordinates
			tensor rot(tensor::zero);
			if (UrelMag != 0)
			{
				vector X(1,0,0);
		        vector R = Urel / UrelMag; 
		        rot = rotationTensor( R, X );
			
				Urel = transform( rot, Urel );
				UTurb = transform( rot, UTurb );
				//Up = transform( rot, Up );
				uiujf = transform( rot, uiujf ); 
			}

			//calculation in transformed coordinates

			//uiujf may not satisfy Cauchy-Scwartz inequality (reduce shear stress to satisfy this property)
			//transformed uiuif may be negative, so we use extra mag() functions to avoid negative sqrt arguments.
			uiujf[0] = max(uiujf[0], uiuiMin_); //max(uiujf[0], ROOTVSMALL); 
			uiujf[3] = max(uiujf[3], uiuiMin_); //max(uiujf[3], ROOTVSMALL);
			uiujf[5] = max(uiujf[5], uiuiMin_); //max(uiujf[5], ROOTVSMALL);
			vector shearStressMax(sqrt(uiujf[0]*uiujf[3]),
								  sqrt(uiujf[0]*uiujf[5]),
								  sqrt(uiujf[3]*uiujf[5]));
			if (mag(uiujf[1])>shearStressMax[0]) uiujf[1]=sign(uiujf[1])*shearStressMax[0];
			if (mag(uiujf[2])>shearStressMax[1]) uiujf[2]=sign(uiujf[2])*shearStressMax[1];
			if (mag(uiujf[4])>shearStressMax[2]) uiujf[4]=sign(uiujf[4])*shearStressMax[2];

			//calculating time sclae tao_s
			vector betai(beta_, 2*beta_, 2*beta_);
			vector unun(uiujf.component(0),uiujf.component(3),uiujf.component(5));
			vector taoL = CL_ * unun / epsilon;

			scalar StE = min(taop/taoL[0] * betai[0], 1000); 

			vector taoT(0, 0, 0);

			for (int i=0; i<3; i++)
			{
				scalar kesiR = min(UrelMag / sqrt(unun[i]), 1e9);
				scalar betaSi = (1-(1-betai[0]*taoL[i]/taoL[0])/(pow(1+StE, 0.4*(1+0.01*StE)) ))* 
					  		    (taoL[0]/taoL[i])*(betai[i]/betai[0]); 
//betaSi = max(min(betaSi, 1e9), 0); max(min(betaSi, 1.0), betai[i]); //max(min(betaSi, 1.0), betai[0]);
				betaSi = max(min(betaSi, 1e9), betai[0]); 
//if (this->owner().db().time().value()>0.0089)  Info << " betaSim = " << betaSi << nl;
//if (btem != betaSi) Info << " betaSi = " << btem << " betaSim = " << betaSi << nl;
//if (this->owner().db().time().value()>0.0017)  Info << "i = " << i << " kesiR = " << kesiR << " betaSi = " << betaSi << " betai[i] = " << betai[i] <<" taoL[i] = " << taoL[i] << " unun[i]=" << unun[i] <<nl;
				scalar tlkS = betaSi/betai[i] * taoL[i]; 
				scalar term = (betaSi*betaSi) * (kesiR*kesiR);
				taoT[i] = tlkS / sqrt(1 + term) + ROOTVSMALL;
			}

			vector ai(vector::zero);  //Eulerian time correlation 
			ai.x() = exp(-dt/taoT.x()); 
			ai.y() = exp(-dt/taoT.y());
			ai.z() = exp(-dt/taoT.z());

			symmTensor sisj; //= symm(tensor::one - ai * ai);
			sisj.xx() = 1 - ai.x() * ai.x();
			sisj.xy() = 1 - ai.x() * ai.y();
			sisj.xz() = 1 - ai.x() * ai.z();
			sisj.yy() = 1 - ai.y() * ai.y();
			sisj.yz() = 1 - ai.y() * ai.z();
			sisj.zz() = 1 - ai.z() * ai.z();

			sisj = cmptMultiply(sisj, uiujf);

			tensor bij = CholskyDecompose(sisj);

			//UTurb = cmptMultiply(ai,UTurb) + bij & randv;
			UTurb = cmptMultiply(ai, UTurb);
			UTurb.x() += bij.xx()*randv.x() + bij.xy()*randv.y() + bij.xz()*randv.z();
			UTurb.y() += bij.yx()*randv.x() + bij.yy()*randv.y() + bij.yz()*randv.z();
			UTurb.z() += bij.zx()*randv.x() + bij.zy()*randv.y() + bij.zz()*randv.z();

			//transforming flactuating velocity back to the global coordinates
			if (UrelMag != 0) UTurb = transform( rot.T(), UTurb );

			//drift correction
			scalar StD = taop / (1./7 * k/epsilon );
			driftCorr  /= (1 + StD);
			UTurb += dt * driftCorr;
			//if (mag(UTurb)>100*mag(Uc)) UTurb = vector::zero;
}

template<class CloudType>
Foam::tensor Foam::CRWDispersionRAS<CloudType>::CholskyDecompose(symmTensor& sisj) 
{
	tensor bij = tensor::zero;

	bij.zz() = sqrt(mag( sisj.zz() ));

	if (bij.zz() > 0) 
	{
		bij.yz() = sisj.yz() / bij.zz();
		bij.xz() = sisj.xz() / bij.zz();
	}

	bij.yy() = sqrt(mag( sisj.yy() - bij.yz()*bij.yz() ));
	if (bij.yy() > 0) bij.xy() = ( sisj.xy() - bij.yz()*bij.xz() ) / bij.yy();

	bij.xx() = sqrt(mag( sisj.xx() - bij.xz()*bij.xz() - bij.xy()*bij.xy() ));

	return bij;
}

template<class CloudType>
void Foam::CRWDispersionRAS<CloudType>::write(Ostream& os) const
{
	DispersionRASModel<CloudType>::write(os);

    os.writeKeyword("ownRf") << ownRf_ << token::END_STATEMENT << endl;
    os.writeKeyword("ownDrift") << ownDrift_ << token::END_STATEMENT
        << endl;
}


// ************************************************************************* //
