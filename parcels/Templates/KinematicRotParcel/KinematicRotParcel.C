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

#include "KinematicRotParcel.H"
#include "torqueSuSp.H"

// * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * * //

template<class ParcelType>
template<class TrackData>
void Foam::KinematicRotParcel<ParcelType>::setCellValues
(
    TrackData& td,
    const scalar dt,
    const label cellI
)
{
    ParcelType::setCellValues(td, dt, cellI);

    //tetIndices tetIs = this->currentTetIndices(); //OmegaC memory

    //Omegac_ = 0.5*td.curlUInterp().interpolate(this->position(), tetIs); //OmegaC memory
}


template<class ParcelType>
template<class TrackData>
void Foam::KinematicRotParcel<ParcelType>::cellValueSourceCorrection
(
    TrackData& td,
    const scalar dt,
    const label cellI
)
{
	ParcelType::cellValueSourceCorrection(td, dt, cellI);
    Omegac_ += td.cloud().OmegaTrans()[cellI]/momentOfInertiaCell(cellI);
}


template<class ParcelType>
template<class TrackData>
void Foam::KinematicRotParcel<ParcelType>::calc
(
    TrackData& td,
    const scalar dt,
    const label cellI
)
{
	ParcelType::calc(td, dt, cellI);

	// Define local properties at beginning of time step
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    const scalar np0 = this->nParticle_;
    const scalar momentOfInertia0 = this->momentOfInertia();

    // Reynolds number
    scalar Re = this->Re(this->U_, this->d_, this->rhoc_, this->muc_); //check

	// Sources
    // ~~~~~~~

	//Explicit momentum source for particle
    vector Somega = vector::zero;

    //Linearised momentum source coefficient
    scalar Spomega = 0.0;

	//Angular Momentum transfer from the particle to the carrier phase
    vector dOmegaTrans = vector::zero;



    // Motion
    // ~~~~~~

    // Calculate new particle angular Momentum
    this->Omega_ = calcOmega(td, dt, cellI, Re, this->muc_, momentOfInertia0, Somega, dOmegaTrans, Spomega);

	//Note2: PIC method should be improves.
    //  Accumulate carrier phase source terms 
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (td.cloud().solution().coupled())
    {
        // Update momentum transfer
        td.cloud().OmegaTrans()[cellI] += np0*dOmegaTrans;

		//Note1: No interaction of omega with carrier phase equations
        // Update momentum transfer coefficient
        //td.cloud().OmegaCoeff()[cellI] += np0*Spomega;
    }
}


//Jafari Added
template<class ParcelType>
template<class TrackData>
const Foam::vector Foam::KinematicRotParcel<ParcelType>::calcOmega
(
    TrackData& td,
    const scalar dt,
    const label cellI,
    const scalar Re,
    const scalar mu,
    const scalar momentOfInertia,
    const vector& Somega,
    vector& dOmegaTrans,
    scalar& Spomega
) //const
{
    typedef typename TrackData::cloudType cloudType;
    typedef typename cloudType::parcelType parcelType;
    typedef typename cloudType::torqueType torqueType;

    const torqueType& torques = td.cloud().torques();

    // Momentum source due to particle forces
    const parcelType& p = static_cast<const parcelType&>(*this);
    scalar taopR(0);
	vector Omegac(vector::zero);
    const torqueSuSp Tcp = torques.calcCoupled(p, dt, momentOfInertia, Re, mu, taopR, Omegac); //taopR and Omegac is calculated
	Omegac_ = Omegac; //OmegaC memory
	//dispersion for Omega td.cloud().dispersionR().tTurbUpdate(tTurbR_, taopR);
    const torqueSuSp Tncp = torques.calcNonCoupled(p, dt, momentOfInertia, Re, mu);
    const torqueSuSp Teff = Tcp + Tncp;
	const scalar momentOfInertiaEff = momentOfInertia; 
	//const scalar momentOfInertiaEff = forces.momentOfInertiaEff(p, momentOfInertia); //for implemeting torques including dwp/dt


    // New particle angular velocity
    //~~~~~~~~~~~~~~~~~~~~~~

    // Update velocity - treat as 3-D
    const vector abp = (Teff.Sp()*Omegac_ + (Teff.Su() + Somega))/momentOfInertiaEff;
    const scalar bp = Teff.Sp()/momentOfInertiaEff;

    Spomega = dt*Teff.Sp();

    IntegrationScheme<vector>::integrationResult Omegares =
        td.cloud().OmegaIntegrator().integrate(Omega_, dt, abp, bp);

    vector Omeganew = Omegares.value();

    // note: Feff.Sp() and Fc.Sp() must be the same
    dOmegaTrans += dt*(Teff.Sp()*(Omegares.average() - Omegac_) - Tcp.Su());

    // Apply correction to velocity and dOmegaTrans for reduced-D cases
    //const polyMesh& mesh = td.cloud().pMesh();
    //meshTools::constrainDirection(mesh, mesh.solutionD(), Omeganew);
    //meshTools::constrainDirection(mesh, mesh.solutionD(), dOmegaTrans);

    return Omeganew;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::KinematicRotParcel<ParcelType>::KinematicRotParcel
(
    const KinematicRotParcel<ParcelType>& p
)
:
    ParcelType(p),
    Omega_(p.Omega_),
	Omegac_(p.Omegac_)
{}


template<class ParcelType>
Foam::KinematicRotParcel<ParcelType>::KinematicRotParcel
(
    const KinematicRotParcel<ParcelType>& p,
    const polyMesh& mesh
)
:
    ParcelType(p, mesh),
    Omega_(p.Omega_),
	Omegac_(p.Omegac_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParcelType>
void Foam::KinematicRotParcel<ParcelType>::transformProperties(const tensor& T)
{
    ParcelType::transformProperties(T);

    Omega_ = transform(T, Omega_);
}


template<class ParcelType>
void Foam::KinematicRotParcel<ParcelType>::transformProperties
(
    const vector& separation
)
{
    ParcelType::transformProperties(separation);
}


// * * * * * * * * * * * * * * IOStream operators  * * * * * * * * * * * * * //

#include "KinematicRotParcelIO.C"

// ************************************************************************* //
