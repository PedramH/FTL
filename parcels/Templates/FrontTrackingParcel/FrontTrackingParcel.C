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

#include "FrontTrackingParcel.H"
//#include "forceSuSp.H"
#include "IntegrationScheme.H"
#include "meshTools.H"

// * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * * //

template<class ParcelType>
template<class TrackData>
void Foam::FrontTrackingParcel<ParcelType>::setCellValues
(
    TrackData& td,
    const scalar dt,
    const label cellI
)
{
    tetIndices tetIs = this->currentTetIndices();

    rhoc_ = td.rhoInterp().interpolate(this->position(), tetIs);

    if (rhoc_ < td.cloud().constProps().rhoMin())
    {
        if (debug)
        {
            WarningIn
            (
                "void Foam::FrontTrackingParcel<ParcelType>::setCellValues"
                "("
                    "TrackData&, "
                    "const scalar, "
                    "const label"
                ")"
            )   << "Limiting observed density in cell " << cellI << " to "
                << td.cloud().constProps().rhoMin() <<  nl << endl;
        }

        rhoc_ = td.cloud().constProps().rhoMin();
    }
  
    Uc_ = td.UInterp().interpolate(this->position(), tetIs);
    this->U_ = Uc_;
/*
        scalar x = this->position().x();
        scalar y = this->position().y();
        scalar z = this->position().z();
        scalar pi = Foam::acos(-1.);
        scalar timeCur = td.cloud().db().time().value();
        this->U_.x() =  2.*Foam::cos(pi*timeCur/8.)*sqr(Foam::sin(pi*x))*Foam::sin(pi*y)*Foam::cos(pi*y);
        this->U_.y() = -2.*Foam::cos(pi*timeCur/8.)*sqr(Foam::sin(pi*y))*Foam::sin(pi*x)*Foam::cos(pi*x);
        this->U_.z() = 0.0;
*/

    muc_ = td.muInterp().interpolate(this->position(), tetIs);
	
    UcMean_ = Uc_;
}


/*
template<class ParcelType>
template<class TrackData>
void Foam::FrontTrackingParcel<ParcelType>::cellValueSourceCorrection
(
    TrackData& td,
    const scalar dt,
    const label cellI
)
{
    Uc_ += td.cloud().UTrans()[cellI]/massCell(cellI);
}
*/

template<class ParcelType>
template<class TrackData>
void Foam::FrontTrackingParcel<ParcelType>::calc
(
    TrackData& td,
    const scalar dt,
    const label cellI
)
{
    // Define local properties at beginning of time step
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    const scalar np0 = nParticle_;
    const scalar mass0 = mass();

    // Reynolds number
    const scalar Re = this->Re(U_, d_, rhoc_, muc_);


    // Sources
    //~~~~~~~~

    // Explicit momentum source for particle
    vector Su = vector::zero;

    // Linearised momentum source coefficient
    scalar Spu = 0.0;

    // Momentum transfer from the particle to the carrier phase
    vector dUTrans = vector::zero;

	vector Up0 = this->U_; //save the value before updating //added ttw
    // Motion
    // ~~~~~~

    // Calculate new particle velocity
    //this->U_ = calcVelocity(td, dt, cellI, Re, muc_, mass0, Su, dUTrans, Spu);

    // *** this is ok
    //tetIndices tetIs = this->currentTetIndices();
    //this->U_ = td.UInterp().interpolate(this->position(), tetIs);
    //        Info << ' ' << "I calculated, U_ and tetIs" << ' ' << td.UInterp().interpolate(this->position(), tetIs)
      //                  << ' ' << this->U_ << endl;

//*** maybe it can be used again
/*
    // Accumulate carrier phase source terms
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (td.cloud().solution().coupled())
    {
        // Update momentum transfer
        td.cloud().UTrans()[cellI] += np0*dUTrans;

        // Update momentum transfer coefficient
        td.cloud().UCoeff()[cellI] += np0*Spu;

		if (td.cloud().solution().turbulenceCoupling()) //added ttw
		{
			symmTensor dRTrans = td.cloud().modulation().update
								 (
									UcMean_,
									Uc_,
									Up0,
									U_,
									dUTrans
								 );

		    // Update Reynold Stress due to solid particles transfer
		    td.cloud().RTrans()[cellI] += np0*dRTrans;
		}
    }
*/

}

//*** maybe it can be used again
/*
template<class ParcelType>
template<class TrackData>
const Foam::vector Foam::FrontTrackingParcel<ParcelType>::calcVelocity
(
    TrackData& td,
    const scalar dt,
    const label cellI,
    const scalar Re,
    const scalar mu,
    const scalar mass,
    const vector& Su,
    vector& dUTrans,
    scalar& Spu
) //const
{
    typedef typename TrackData::cloudType cloudType;
    typedef typename cloudType::parcelType parcelType;
    typedef typename cloudType::forceType forceType;

    const forceType& forces = td.cloud().forces();

    // Momentum source due to particle forces
    const parcelType& p = static_cast<const parcelType&>(*this);
	//const forceSuSp Fcp = forces.calcCoupled(p, dt, mass, Re, mu);

	scalar taop(0);
    const forceSuSp Fcp = forces.calcCoupled(p, dt, mass, Re, mu, taop); //taop is calculated
	td.cloud().dispersion().tTurbUpdate(tTurb_, taop);
	//Info << "tTurb = " << tTurb_ << "taop = " << taop << nl;

    const forceSuSp Fncp = forces.calcNonCoupled(p, dt, mass, Re, mu);
    const forceSuSp Feff = Fcp + Fncp;
    const scalar massEff = forces.massEff(p, mass);


    // New particle velocity
    //~~~~~~~~~~~~~~~~~~~~~~

    // Update velocity - treat as 3-D
    const vector abp = (Feff.Sp()*Uc_ + (Feff.Su() + Su))/massEff;
    const scalar bp = Feff.Sp()/massEff;

    Spu = dt*Feff.Sp();

    IntegrationScheme<vector>::integrationResult Ures =
        td.cloud().UIntegrator().integrate(U_, dt, abp, bp);

    vector Unew = Ures.value();

    // note: Feff.Sp() and Fc.Sp() must be the same
    dUTrans += dt*(Feff.Sp()*(Ures.average() - Uc_) - Fcp.Su());

    // Apply correction to velocity and dUTrans for reduced-D cases
    const polyMesh& mesh = td.cloud().pMesh();
    meshTools::constrainDirection(mesh, mesh.solutionD(), Unew);
    meshTools::constrainDirection(mesh, mesh.solutionD(), dUTrans);

    return Unew;
}
*/

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::FrontTrackingParcel<ParcelType>::FrontTrackingParcel
(
    const FrontTrackingParcel<ParcelType>& p
)
:
    ParcelType(p),
    active_(p.active_),
    typeId_(p.typeId_),
    nParticle_(p.nParticle_),
    d_(p.d_),
    dTarget_(p.dTarget_),
    U_(p.U_),
    rho_(p.rho_),
    age_(p.age_),
    tTurb_(p.tTurb_),
    UTurb_(p.UTurb_),
    currentIndex_(p.currentIndex_),//Jafari added
    bubbleIndex_(p.bubbleIndex_),//Jafari added
    rhoc_(p.rhoc_),
    Uc_(p.Uc_),
    muc_(p.muc_)
{}


template<class ParcelType>
Foam::FrontTrackingParcel<ParcelType>::FrontTrackingParcel
(
    const FrontTrackingParcel<ParcelType>& p,
    const polyMesh& mesh
)
:
    ParcelType(p, mesh),
    active_(p.active_),
    typeId_(p.typeId_),
    nParticle_(p.nParticle_),
    d_(p.d_),
    dTarget_(p.dTarget_),
    U_(p.U_),
    rho_(p.rho_),
    age_(p.age_),
    tTurb_(p.tTurb_),
    UTurb_(p.UTurb_),
    currentIndex_(p.currentIndex_),//Jafari added
    bubbleIndex_(p.bubbleIndex_),//Jafari added
    rhoc_(p.rhoc_),
    Uc_(p.Uc_),
    muc_(p.muc_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParcelType>
template<class TrackData>
bool Foam::FrontTrackingParcel<ParcelType>::move
(
    TrackData& td,
    const scalar trackTime
)
{
    typename TrackData::cloudType::parcelType& p =
        static_cast<typename TrackData::cloudType::parcelType&>(*this);

    td.switchProcessor = false;
    td.keepParticle = true;

    const polyMesh& mesh = td.cloud().pMesh();
    const polyBoundaryMesh& pbMesh = mesh.boundaryMesh();
    const scalarField& cellLengthScale = td.cloud().cellLengthScale();
    const scalar maxCo = td.cloud().solution().maxCo();

    scalar tEnd = (1.0 - p.stepFraction())*trackTime;
    const scalar dtMax = tEnd;

    bool moving = true;
/*
            Info << ' ' << "p.stepFraction()" << ' ' << p.stepFraction() << endl;
            Info << ' ' << "trackTime" << ' ' << trackTime << endl;

            Info << ' ' << "tEnd" << ' ' << tEnd << endl;
            Info << ' ' << "td.keepParticle" << ' ' << td.keepParticle << endl;
            Info << ' ' << "td.switchProcessor" << ' ' << td.switchProcessor << endl;
*/

    while (td.keepParticle && !td.switchProcessor && tEnd > ROOTVSMALL)
    {

            //Info << ' ' << "I'm in while .....loop" << endl;
        // Apply correction to position for reduced-D cases
        meshTools::constrainToMeshCentre(mesh, p.position());

        const point start(p.position());

        // Set the Lagrangian time-step
        scalar dt = min(dtMax, tEnd);

        // Cache the parcel current cell as this will change if a face is hit
        const label cellI = p.cell();

        //Jafrai added for uodating the front point varaiables like velocity, pressure, temprature,...
        this->U_ =td.UInterp().interpolate( this->position(), this->currentTetIndices());
/*
        scalar x = p.position().x();
        scalar y = p.position().y();
        scalar z = p.position().z();
        scalar pi = Foam::acos(-1.);
        scalar timeCur = td.cloud().db().time().value();
        this->U_.x() =  2.*Foam::cos(pi*timeCur/8.)*sqr(Foam::sin(pi*x))*Foam::sin(pi*y)*Foam::cos(pi*y);
        this->U_.y() = -2.*Foam::cos(pi*timeCur/8.)*sqr(Foam::sin(pi*y))*Foam::sin(pi*x)*Foam::cos(pi*x);
        this->U_.z() = 0.0;
*/
        const scalar magU = mag(U_);
        //Info << ' ' << "My p.active()/moving/magU are " << p.active() << ' ' << moving << ' ' << magU << endl;

        if (p.active() && moving && (magU > ROOTVSMALL))
        {
            const scalar d = dt*magU;
            const scalar dCorr = min(d, maxCo*cellLengthScale[cellI]);
            dt *=
                dCorr/d
               *p.trackToFace(p.position() + dCorr*U_/magU, td);
            //Info << ' ' << "am I here in  transferring loop with delX: " << (this->position() -this->shadowPos()) << endl;
        }

        tEnd -= dt;

        scalar newStepFraction = 1.0 - tEnd/trackTime;

        if
        (
            mag(p.stepFraction() - newStepFraction)
          < particle::minStepFractionTol
        )
        {
            moving = false;
        }

        p.stepFraction() = newStepFraction;

        // Avoid problems with extremely small timesteps
        if (dt > ROOTVSMALL)
        {
            // Update cell based properties
            //p.setCellValues(td, dt, cellI);

            //if (td.cloud().solution().cellValueSourceCorrection())
            //{
                //p.cellValueSourceCorrection(td, dt, cellI);
            //}
            //Info << ' ' << "I'm before calc" << endl;
            //p.calc(td, dt, cellI);
        }

        if (p.onBoundary() && td.keepParticle)
        {
            if (isA<processorPolyPatch>(pbMesh[p.patch(p.face())]))
            {
                td.switchProcessor = true;
            }
        }

        p.age() += dt;

        //td.cloud().functions().postMove(p, cellI, dt, start, td.keepParticle);
    }

    return td.keepParticle;
}


template<class ParcelType>
template<class TrackData>
void Foam::FrontTrackingParcel<ParcelType>::hitFace(TrackData& td)
{
    typename TrackData::cloudType::parcelType& p =
        static_cast<typename TrackData::cloudType::parcelType&>(*this);

    //td.cloud().functions().postFace(p, p.face(), td.keepParticle);
}


template<class ParcelType>
void Foam::FrontTrackingParcel<ParcelType>::hitFace(int& td)
{}


template<class ParcelType>
template<class TrackData>
bool Foam::FrontTrackingParcel<ParcelType>::hitPatch
(
    const polyPatch& pp,
    TrackData& td,
    const label patchI,
    const scalar trackFraction,
    const tetIndices& tetIs
)
{
    typename TrackData::cloudType::parcelType& p =
        static_cast<typename TrackData::cloudType::parcelType&>(*this);

    // Invoke post-processing model
/*
    td.cloud().functions().postPatch
    (
        p,
        pp,
        trackFraction,
        tetIs,
        td.keepParticle
    );
*/

/*
    // Invoke surface film model
    if (td.cloud().surfaceFilm().transferParcel(p, pp, td.keepParticle))
    {
        // All interactions done
        return true;
    }
    else
    {
        // Invoke patch interaction model
        return td.cloud().patchInteraction().correct
        (
            p,
            pp,
            td.keepParticle,
            trackFraction,
            tetIs
        );
    }
*/
        return false;
}


template<class ParcelType>
template<class TrackData>
void Foam::FrontTrackingParcel<ParcelType>::hitProcessorPatch
(
    const processorPolyPatch&,
    TrackData& td
)
{
    td.switchProcessor = true;
}


template<class ParcelType>
template<class TrackData>
void Foam::FrontTrackingParcel<ParcelType>::hitWallPatch
(
    const wallPolyPatch& wpp,
    TrackData& td,
    const tetIndices&
)
{
    // Wall interactions handled by generic hitPatch function
}


template<class ParcelType>
template<class TrackData>
void Foam::FrontTrackingParcel<ParcelType>::hitPatch
(
    const polyPatch&,
    TrackData& td
)
{
    td.keepParticle = false;
}


template<class ParcelType>
void Foam::FrontTrackingParcel<ParcelType>::transformProperties(const tensor& T)
{
    ParcelType::transformProperties(T);

    U_ = transform(T, U_);
}


template<class ParcelType>
void Foam::FrontTrackingParcel<ParcelType>::transformProperties
(
    const vector& separation
)
{
    ParcelType::transformProperties(separation);
}


template<class ParcelType>
Foam::scalar Foam::FrontTrackingParcel<ParcelType>::wallImpactDistance
(
    const vector&
) const
{
    return 0.5*d_;
}


// * * * * * * * * * * * * * * IOStream operators  * * * * * * * * * * * * * //

#include "FrontTrackingParcelIO.C"

// ************************************************************************* //
