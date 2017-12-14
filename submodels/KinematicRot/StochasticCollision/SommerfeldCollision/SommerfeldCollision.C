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

#include "SommerfeldCollision.H"
#include "mathematicalConstants.H"
#include "turbulenceModel.H" //changed
//#include "SLGThermo.H"


using namespace Foam::constant;

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class CloudType>
void Foam::SommerfeldCollision<CloudType>::cacheFields(const bool store)
{
    if (store)
    {
        const fvMesh& mesh = this->owner().mesh();
        const word& cloudName = this->owner().name();

        const DimensionedField<vector, volMesh>& UpMean =
            mesh.lookupObject<DimensionedField<vector, volMesh> >
            (
                cloudName + ":UPStat"
            );
		UpMean_ = &UpMean;
		const DimensionedField<vector, volMesh>& upup =
            mesh.lookupObject<DimensionedField<vector, volMesh> >
            (
                cloudName + ":uPuPStat"
            );
		upup_ = &upup;
		const DimensionedField<scalar, volMesh>& NpMean =
            mesh.lookupObject<DimensionedField<scalar, volMesh> >
            (
                cloudName + ":NPStat"
            );
		NpMean_ = &NpMean;
		const DimensionedField<scalar, volMesh>& dpMean =
            mesh.lookupObject<DimensionedField<scalar, volMesh> >
            (
                cloudName + ":DPStat"
            );
		dpMean_ = &dpMean;
		const DimensionedField<scalar, volMesh>& dpdp =
            mesh.lookupObject<DimensionedField<scalar, volMesh> >
            (
                cloudName + ":dPdPStat"
            );
		dpdp_ = &dpdp;
		const DimensionedField<vector, volMesh>& OmegapMean =
            mesh.lookupObject<DimensionedField<vector, volMesh> >
            (
                cloudName + ":OmegaPStat"
            );
		OmegapMean_ = &OmegapMean;
    }
    else
    {
        UpMean_ = NULL;
        upup_ = NULL;
		NpMean_ = NULL;
		dpMean_ = NULL;
		dpdp_ = NULL;
		OmegapMean_ = NULL;
    }
	//from DispersionRASModel
	if (store)
    {
        tmp<volScalarField> tk = this->kModel();
        if (tk.isTmp())
        {
            kPtr_ = tk.ptr();
            ownK_ = true;
        }
        else
        {
            kPtr_ = tk.operator->();
            ownK_ = false;
        }

        tmp<volScalarField> tepsilon = this->epsilonModel();
        if (tepsilon.isTmp())
        {
            epsilonPtr_ = tepsilon.ptr();
            ownEpsilon_ = true;
        }
        else
        {
            epsilonPtr_ = tepsilon.operator->();
            ownEpsilon_ = false;
        }

		tmp<volScalarField> trho = this->rhoModel();
        if (trho.isTmp())
        {
            rhoPtr_ = trho.ptr();
            ownRho_ = true;
        }
        else
        {
            rhoPtr_ = trho.operator->();
            ownRho_ = false;
        }

		tmp<volScalarField> tmu = this->muModel();
        if (tmu.isTmp())
        {
            muPtr_ = tmu.ptr();
            ownMu_ = true;
        }
        else
        {
            muPtr_ = tmu.operator->();
            ownMu_ = false;
        }
    }
    else
    {
        if (ownK_ && kPtr_)
        {
            deleteDemandDrivenData(kPtr_);
            ownK_ = false;
        }
        if (ownEpsilon_ && epsilonPtr_)
        {
            deleteDemandDrivenData(epsilonPtr_);
            ownEpsilon_ = false;
        }
		if (ownRho_ && rhoPtr_)
        {
            deleteDemandDrivenData(rhoPtr_);
            ownRho_ = false;
        }
		if (ownMu_ && muPtr_)
        {
            deleteDemandDrivenData(muPtr_);
            ownMu_ = false;
        }
    }
}


template<class CloudType>
Foam::tmp<Foam::volScalarField>
Foam::SommerfeldCollision<CloudType>::kModel() const
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
        return model.k();
    }
    else
    {
        FatalErrorIn
        (
            "Foam::tmp<Foam::volScalarField>"
            "Foam::DispersionRASModel<CloudType>::kModel() const"
        )
            << "Turbulence model not found in mesh database" << nl
            << "Database objects include: " << obr.sortedToc()
            << abort(FatalError);

        return tmp<volScalarField>(NULL);
    }
}


template<class CloudType>
Foam::tmp<Foam::volScalarField>
Foam::SommerfeldCollision<CloudType>::epsilonModel() const
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
        return model.epsilon();
    }
    else
    {
        FatalErrorIn
        (
            "Foam::tmp<Foam::volScalarField>"
            "Foam::DispersionRASModel<CloudType>::epsilonModel() const"
        )
            << "Turbulence model not found in mesh database" << nl
            << "Database objects include: " << obr.sortedToc()
            << abort(FatalError);

        return tmp<volScalarField>(NULL);
    }
}


template<class CloudType>
Foam::tmp<Foam::volScalarField>
Foam::SommerfeldCollision<CloudType>::rhoModel() const
{
    const objectRegistry& obr = this->owner().mesh();

    if (obr.foundObject<volScalarField>("rho"))
    {
		const volScalarField& rhoc =
            obr.lookupObject<volScalarField>("rho");

        return rhoc;
    }
    else
    {
        FatalErrorIn
        (
            "Foam::tmp<Foam::volScalarField>"
            "Foam::SommerfeldCollision<CloudType>::rhoModel() const"
        )
            << "carrier phase fields not found in mesh database" << nl
            << "Database objects include: " << obr.sortedToc()
            << abort(FatalError);

        return tmp<volScalarField>(NULL);
    }
}


template<class CloudType>
Foam::tmp<Foam::volScalarField>
Foam::SommerfeldCollision<CloudType>::muModel() const
{
    const objectRegistry& obr = this->owner().mesh();

    if (obr.foundObject<volScalarField>("mu"))
    {
		const volScalarField& muc =
            obr.lookupObject<volScalarField>("mu");

        return muc;
    }
    else
    {
        FatalErrorIn
        (
            "Foam::tmp<Foam::volScalarField>"
            "Foam::SommerfeldCollision<CloudType>::muModel() const"
        )
            << "carrier phase fields not found in mesh database" << nl
            << "Database objects include: " << obr.sortedToc()
            << abort(FatalError);

        return tmp<volScalarField>(NULL);
    }
}


template<class CloudType>
void Foam::SommerfeldCollision<CloudType>::collide(const scalar dt)
{
	cacheFields(true);

	int numberOfCollisions = 0;
    forAllIter(typename CloudType, this->owner(), iter)
    {
		parcelType& p = iter();
		bool doCollide = collideParcels(dt, p);
		if (doCollide) numberOfCollisions++;
    }
	Info << " " << numberOfCollisions << " binary collisions occur between parcels." << nl;

	cacheFields(false);
}


template<class CloudType>
bool Foam::SommerfeldCollision<CloudType>::collideParcels
(
    const scalar dt,
    parcelType& p
)
{
    const label cellI = p.cell();

	const scalar d1 = p.d();

	const scalar d2Mean = dpMean_->field()[cellI];
	const scalar d2rms = sqrt(mag(dpdp_->field()[cellI]));
	scalar RN = ranGen_.scalar01();
	const scalar d2 = RN*d2rms + d2Mean;

	const scalar m1 = p.mass();
	const scalar rho2 = this->owner().constProps().rho0(); //p.rho();
	const scalar m2 = rho2*mathematical::pi/6*pow(d2,3); 

    if ((m1 < ROOTVSMALL) || (m2 < ROOTVSMALL)) return false;
	
	vector& U1 = p.U(); //changed by pointer
	const vector UpMean = UpMean_->field()[cellI];
	const vector UpStd = upup_->field()[cellI];
	const vector Uprms(sqrt(mag(UpStd.x())), sqrt(mag(UpStd.y())), sqrt(mag(UpStd.z())));   //sqrt(upup_->field()[cellI]);

	const scalar muc = muPtr_->internalField()[cellI];
	scalar taopSt = rho2*sqr(d2)/(18*muc + ROOTVSMALL);
	const scalar k = kPtr_->internalField()[cellI] + ROOTVSMALL;
    const scalar epsilon = epsilonPtr_->internalField()[cellI] + ROOTVSMALL;
	scalar TE = 0.3*k/epsilon;
	scalar RSt = exp(-0.55*pow(taopSt/TE,0.4));
	vector randv;
	randv.x() = ranGen_.GaussNormal();
	randv.y() = ranGen_.GaussNormal();
	randv.z() = ranGen_.GaussNormal();
	vector U2 = UpMean + (RSt*(U1-UpMean) + sqrt(1-RSt*RSt)*cmptMultiply(Uprms,randv));

	vector G0 = U1 - U2;
	scalar magG0 = mag(G0);
	const scalar np = NpMean_->field()[cellI];
	scalar fc = mathematical::pi*pow((d1+d2)/2,2)*magG0*np;
	scalar Pcoll = 1-exp(-fc*dt);
 
	RN = ranGen_.scalar01();
	if (RN >= Pcoll) return false; //it is better to limit the time step by C*toac toac=1/fc, C<<1 (C=0.1)

	//position of fictious particel, 2.
	scalar L;
	do 
	{
		scalar XX = ranGen_.scalar01();
		scalar YY = ranGen_.scalar01();
		L = sqrt(XX*XX+YY*YY);
	} while (L >= 1);

	if (repulsiveCorrection_) 
	{
		scalar La = L * (d1+d2)/2;
		scalar dMax = max(d1, d2);
		scalar dMin = min(d1, d2);
		if (dMax/dMin > 10)
		{
			const scalar rhoc = rhoPtr_->internalField()[cellI];
			scalar Yc;
			if (dMax == d1)
			{
				Yc = limitingDistance(muc,rhoc,d1,U1,rho2,d2,U2);
			}
			else
			{
				Yc = limitingDistance(muc,rhoc,d2,U2,p.rho(),d1,U1);
			}
			if (La >= Yc) return false; //condition4
		}
	}

	//calculating collision normal vector from 1 to 2

	scalar phi = asin(L);
	RN = ranGen_.scalar01();
	scalar psi = 2*mathematical::pi*RN;

	vector xr21(sin(phi)*cos(psi), sin(phi)*sin(psi), cos(phi)); // xr21=x2-x1 in relative coordinates
	xr21 *= (d1+d2)/2;
	//transforming relative position back to the global coordinates
	tensor rot(tensor::zero);
	vector X(0,0,1);
    vector R = G0 / magG0;  //magG0 > 0
    rot = rotationTensor( R, X );
	vector xrCol = -transform( rot.T(), xr21 ); //xrCol=x1-x2

	vector nC = -xrCol/(mag(xrCol) + ROOTVSMALL);

	//calculating collision outcome
	
	vector& Omega1 = p.Omega(); //changed by pointer
	vector Omega2 = OmegapMean_->field()[cellI];

	vector Gc0 = G0 + 0.5*d1*(Omega1 ^ nC) + 0.5*d2*(Omega2 ^ nC);
    vector Gct0 = Gc0-(Gc0 & nC)*nC;
	scalar magGct0 = mag(Gct0);
    vector tC = Gct0/(magGct0 + ROOTVSMALL);

	if ( (nC & G0)/(mag(Gct0) + ROOTVSMALL) < 2.0/( 7.0*mu_*(1+e_)) ) // sliding collision 
    {
        scalar coe = (nC & G0)*(1+e_);
		scalar coe1 = coe*m2/(m1+m2);
		//scalar coe2 = coe*m1/(m1+m2);
		U1 = U1 - (nC+mu_*tC)*coe1;
		//U2 = U2 + (nC+mu_*tC)*coe2;
		Omega1 = Omega1 - (5.0/d1)*coe1*(nC ^ tC)*mu_;
		//Omega2 = Omega2 - (5.0/d2)*coe2*(nC ^ tC)*mu_;
    }
    else 
    {
		vector vc = (1+e_)*(nC & G0)*nC + (2.0/7.0)*mag(Gct0)*tC;
		U1= U1 - vc*m2/(m1+m2);
		//U2= U2 + vc*m1/(m1+m2);
		Omega1= Omega1 - (10.0/7.0)*magGct0/d1*(nC ^ tC)*m2/(m1+m2);
		//Omega2= Omega2 - (10.0/7.0)*magGct0/d2*(nC ^ tC)*m1/(m1+m2);
    }

    return true;
}


template<class CloudType>
Foam::scalar Foam::SommerfeldCollision<CloudType>::limitingDistance 
(
    const scalar muc,
    const scalar rhoc,
	const scalar dL,
	const vector UL,
	const scalar rhoS,
	const scalar dS,
	const vector US
)
{
	scalar Str = rhoS*sqr(dS)*mag(US-UL)/(18*muc*dL + ROOTVSMALL);
	scalar ReL = mag(UL)*dL*rhoc/(muc + ROOTVSMALL);
	scalar a, b;
	if (ReL <= 1)
	{
		a = 0.65;
		b = 3.7;
	}
	else if (ReL < 30)
	{
		a = 1.24;
		b = 1.95;
	}
	else if (ReL < 50)
	{
		a = 1.03;
		b = 2.07;
	}
	else if (ReL < 80)
	{
		a = 0.506;
		b = 1.84;
	}
	else
	{
		a = 0.25;
		b = 2.0;
	}
	return dL/2*pow(Str/(Str+a), b/2);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::SommerfeldCollision<CloudType>::SommerfeldCollision
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName
)
:
    StochasticCollisionModel<CloudType>(dict, owner, modelName),
	ranGen_(label(0)),
	e_(readScalar(this->coeffDict().lookup("e"))),
    mu_(readScalar(this->coeffDict().lookup("mu"))),
	repulsiveCorrection_(this->coeffDict().lookupOrDefault("repulsiveCorrection",false)),
	UpMean_(NULL),
    upup_(NULL),
	NpMean_(NULL),
	dpMean_(NULL),
	dpdp_(NULL),
	OmegapMean_(NULL),
	kPtr_(NULL),
    ownK_(false),
    epsilonPtr_(NULL),
    ownEpsilon_(false),
	rhoPtr_(NULL),
    ownRho_(false),
    muPtr_(NULL),
    ownMu_(false)
    /*liquids_
    (
        owner.db().template lookupObject<SLGThermo>("SLGThermo").liquids()
    ),
    coalescence_(this->coeffDict().lookup("coalescence"))*/
{
	//change cloud flag to compute statistics group 1 and 2
	owner.preColStat1Comp() = true;
	owner.preColStat2Comp() = true;
}


template<class CloudType>
Foam::SommerfeldCollision<CloudType>::SommerfeldCollision
(
    SommerfeldCollision<CloudType>& cm
)
:
    StochasticCollisionModel<CloudType>(cm),
	ranGen_(cm.ranGen_),
	e_(cm.e_),
    mu_(cm.mu_),
	repulsiveCorrection_(cm.repulsiveCorrection_),
	UpMean_(cm.UpMean_),
    upup_(cm.upup_),
	NpMean_(cm.NpMean_),
	dpMean_(cm.dpMean_),
	dpdp_(cm.dpdp_),
	OmegapMean_(cm.OmegapMean_),
	kPtr_(cm.kPtr_),
    ownK_(cm.ownK_),
    epsilonPtr_(cm.epsilonPtr_),
    ownEpsilon_(cm.ownEpsilon_),
	rhoPtr_(cm.rhoPtr_),
    ownRho_(cm.ownRho_),
    muPtr_(cm.muPtr_),
    ownMu_(cm.ownMu_)
    /*liquids_(cm.liquids_),
    coalescence_(cm.coalescence_)*/
{
	cm.ownK_ = false;
    cm.ownEpsilon_ = false;
	cm.ownRho_ = false;
    cm.ownMu_ = false;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::SommerfeldCollision<CloudType>::~SommerfeldCollision()
{}


// ************************************************************************* //
