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

#include "WDeterministicCollision.H"
//#include "SLGThermo.H"


using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //


template<class CloudType>
void Foam::WDeterministicCollision<CloudType>::cacheFields(const bool store)
{
	if (repulsiveCorrection_) 
	{
		//from DispersionRASModel
		if (store)
		{
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
}


template<class CloudType>
Foam::tmp<Foam::volScalarField>
Foam::WDeterministicCollision<CloudType>::rhoModel() const
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
            "Foam::WDeterministicCollision<CloudType>::rhoModel() const"
        )
            << "carrier phase fields not found in mesh database" << nl
            << "Database objects include: " << obr.sortedToc()
            << abort(FatalError);

        return tmp<volScalarField>(NULL);
    }
}


template<class CloudType>
Foam::tmp<Foam::volScalarField>
Foam::WDeterministicCollision<CloudType>::muModel() const
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
            "Foam::WDeterministicCollision<CloudType>::muModel() const"
        )
            << "carrier phase fields not found in mesh database" << nl
            << "Database objects include: " << obr.sortedToc()
            << abort(FatalError);

        return tmp<volScalarField>(NULL);
    }
}


template<class CloudType>
void Foam::WDeterministicCollision<CloudType>::collide(const scalar dt)
{
	cacheFields(true);

	int numberOfCollisions = 0;
    label i = 0;
    forAllIter(typename CloudType, this->owner(), iter1)
    {
        label j = 0;
        forAllIter(typename CloudType, this->owner(), iter2)
        {
            if (j > i)
            {
                parcelType& p1 = iter1();
                parcelType& p2 = iter2();

				scalar m1 = p1.mass(); //p1.nParticle()*p1.mass();
                scalar m2 = p2.mass(); //p2.nParticle()*p2.mass();

				bool doCollide = collideParcels(dt, p1, p2, m1, m2);
				if (doCollide) numberOfCollisions++;
            }
            j++;
        }
        i++;
    }
	Info << " " << numberOfCollisions << " binary collisions occur between parcels." << nl;

	cacheFields(false);
}


template<class CloudType>
bool Foam::WDeterministicCollision<CloudType>::collideParcels
(
    const scalar dt,
    parcelType& p1,
    parcelType& p2,
    scalar& m1,
    scalar& m2
)
{
    const label cell1 = p1.cell();
    const label cell2 = p2.cell();

	labelList neighbours1 = this->owner().mesh().cellCells()[cell1];
	//for adding corner cells 
	labelList corners1;
	forAll(neighbours1,j)
	{
		corners1.append(this->owner().mesh().cellCells()[neighbours1[j]]);
	}
	neighbours1.append(corners1); //note: there are repeatitive cells in neighbours1, but it does not matter.	

    // check if parcels belong to same or neighbor cells
    if ((m1 < ROOTVSMALL) || (m2 < ROOTVSMALL))
    {
		return false;
	}
	else if (cell1 != cell2)
	{
		bool flag = true;
		forAll(neighbours1,j)
		{
			if (neighbours1[j] == cell2) 
			{
				flag = false;
				break;
			}
		}
		if (flag) return false;
    }

	vector U1 = p1.U();
	vector U2 = p2.U();
	vector G0 = U1 - U2;
	vector xr0 = (p1.position()-p2.position()) - dt * G0;
	scalar vAlign = xr0 & G0;

	if (vAlign >= 0) return false; //condition1
	
	scalar magG0 = mag(G0);
	scalar dtMin = -vAlign / sqr(magG0); //mag(G0) > 0
 
	if (dtMin > dt) return false; //condition2

	const scalar d1 = p1.d();
    const scalar d2 = p2.d();
	scalar nP1 = p1.nParticle();
    scalar nP2 = p2.nParticle();
	scalar nMax = max(nP1, nP2);
	scalar sigma = sqr((d1+d2)/2) * nMax; //parcel

	scalar f = (magSqr(xr0)-sigma)/sqr(dtMin*magG0); //dtMin > 0

	if (f >= 1) return false; //condition3

	if (repulsiveCorrection_) 
	{
		scalar dMax = max(d1, d2);
		scalar dMin = min(d1, d2);
		if (dMax/dMin > 10)
		{
			scalar La = mag(xr0 + dtMin * G0)/sqrt(nMax);
			scalar Yc;
			if (dMax == d1)
			{
				const scalar rhoc = rhoPtr_->internalField()[cell1];
				const scalar muc = muPtr_->internalField()[cell1];
				Yc = limitingDistance(muc,rhoc,d1,U1,p2.rho(),d2,U2);
			}
			else
			{
				const scalar rhoc = rhoPtr_->internalField()[cell2];
				const scalar muc = muPtr_->internalField()[cell2];
				Yc = limitingDistance(muc,rhoc,d2,U2,p1.rho(),d1,U1);
			}
			if (La >= Yc) return false; //condition4
		}
	}

	scalar dtCol = dtMin * (1-sqrt(1-f));
	vector xrCol = xr0 + dtCol * G0;

	vector nC = -xrCol/(mag(xrCol) + ROOTVSMALL);
	
	vector Omega1 = p1.Omega();
    vector Omega2 = p2.Omega();

	vector Gc0 = G0 + 0.5*d1*(Omega1 ^ nC) + 0.5*d2*(Omega2 ^ nC);
    vector Gct0 = Gc0-(Gc0 & nC)*nC;
	scalar magGct0 = mag(Gct0);
    vector tC = Gct0/(magGct0 + ROOTVSMALL);

	if ( (nC & G0)/(mag(Gct0) + ROOTVSMALL) < 2.0/( 7.0*mu_*(1+e_)) ) // sliding collision 
    {
        scalar coe = (nC & G0)*(1+e_);
		scalar coe1 = coe*m2/(m1+m2);
		scalar coe2 = coe*m1/(m1+m2);
		U1 = U1 - (nC+mu_*tC)*coe1;
		U2 = U2 + (nC+mu_*tC)*coe2;
		Omega1 = Omega1 - (5.0/d1)*coe1*(nC ^ tC)*mu_;
		Omega2 = Omega2 - (5.0/d2)*coe2*(nC ^ tC)*mu_;
    }
    else 
    {
		vector vc = (1+e_)*(nC & G0)*nC + (2.0/7.0)*mag(Gct0)*tC;
		U1= U1 - vc*m2/(m1+m2);
		U2= U2 + vc*m1/(m1+m2);
		Omega1= Omega1 - (10.0/7.0)*magGct0/d1*(nC ^ tC)*m2/(m1+m2);
		Omega2= Omega2 - (10.0/7.0)*magGct0/d2*(nC ^ tC)*m1/(m1+m2);
    }

	if (nP1 < nP2) //nMin collisions occur
	{
        p1.U() = U1;
        p2.U() = (nP1*U2 + (nP2 - nP1)*p2.U())/nP2;
		p1.Omega() = Omega1;
        p2.Omega() = (nP1*Omega2 + (nP2 - nP1)*p2.Omega())/nP2;
    }
    else
    {
        p1.U() = (nP2*U1 + (nP1 - nP2)*p1.U())/nP1;
        p2.U() = U2;
		p1.Omega() = (nP2*Omega1 + (nP1 - nP2)*p1.Omega())/nP1;
        p2.Omega() = Omega2;
    }

    return true;
}


template<class CloudType>
Foam::scalar Foam::WDeterministicCollision<CloudType>::limitingDistance 
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
Foam::WDeterministicCollision<CloudType>::WDeterministicCollision
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName
)
:
    StochasticCollisionModel<CloudType>(dict, owner, modelName),
	e_(readScalar(this->coeffDict().lookup("e"))),
    mu_(readScalar(this->coeffDict().lookup("mu"))),
	repulsiveCorrection_(this->coeffDict().lookupOrDefault("repulsiveCorrection",false)),
	rhoPtr_(NULL),
    ownRho_(false),
    muPtr_(NULL),
    ownMu_(false)
    /*liquids_
    (
        owner.db().template lookupObject<SLGThermo>("SLGThermo").liquids()
    ),
    coalescence_(this->coeffDict().lookup("coalescence"))*/
{}


template<class CloudType>
Foam::WDeterministicCollision<CloudType>::WDeterministicCollision
(
    WDeterministicCollision<CloudType>& cm
)
:
    StochasticCollisionModel<CloudType>(cm),
	e_(cm.e_),
    mu_(cm.mu_),
	repulsiveCorrection_(cm.repulsiveCorrection_),
	rhoPtr_(cm.rhoPtr_),
    ownRho_(cm.ownRho_),
    muPtr_(cm.muPtr_),
    ownMu_(cm.ownMu_)
    /*liquids_(cm.liquids_),
    coalescence_(cm.coalescence_)*/
{
	cm.ownRho_ = false;
    cm.ownMu_ = false;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::WDeterministicCollision<CloudType>::~WDeterministicCollision()
{}


// ************************************************************************* //
