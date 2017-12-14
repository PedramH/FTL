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

#include "fvmSup.H"
#include "SortableList.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
inline const Foam::FrontTrackingCloud<CloudType>&
Foam::FrontTrackingCloud<CloudType>::cloudCopy() const
{
    return cloudCopyPtr_();
}


template<class CloudType>
inline const Foam::fvMesh& Foam::FrontTrackingCloud<CloudType>::mesh() const
{
    return mesh_;
}


template<class CloudType>
inline const Foam::IOdictionary&
Foam::FrontTrackingCloud<CloudType>::particleProperties() const
{
    return particleProperties_;
}


template<class CloudType>
inline const Foam::IOdictionary&
Foam::FrontTrackingCloud<CloudType>::outputProperties() const
{
    return outputProperties_;
}


template<class CloudType>
inline Foam::IOdictionary& Foam::FrontTrackingCloud<CloudType>::outputProperties()
{
    return outputProperties_;
}


template<class CloudType>
inline const Foam::cloudSolution&
Foam::FrontTrackingCloud<CloudType>::solution() const
{
    return solution_;
}


template<class CloudType>
inline Foam::cloudSolution& Foam::FrontTrackingCloud<CloudType>::solution()
{
    return solution_;
}


template<class CloudType>
inline const typename CloudType::particleType::constantProperties&
Foam::FrontTrackingCloud<CloudType>::constProps() const
{
    return constProps_;
}


template<class CloudType>
inline typename CloudType::particleType::constantProperties&
Foam::FrontTrackingCloud<CloudType>::constProps()
{
    return constProps_;
}


template<class CloudType>
inline const Foam::dictionary&
Foam::FrontTrackingCloud<CloudType>::subModelProperties() const
{
    return subModelProperties_;
}


template<class CloudType>
inline const Foam::volScalarField& Foam::FrontTrackingCloud<CloudType>::rho() const
{
    return rho_;
}


template<class CloudType>
inline const Foam::volVectorField& Foam::FrontTrackingCloud<CloudType>::U() const
{
    return U_;
}


template<class CloudType>
inline const Foam::volScalarField& Foam::FrontTrackingCloud<CloudType>::mu() const
{
    return mu_;
}

// Jafari added
template<class CloudType>
inline const Foam::volScalarField& Foam::FrontTrackingCloud<CloudType>::p() const
{
    return p_;
}


template<class CloudType>
inline const Foam::dimensionedVector& Foam::FrontTrackingCloud<CloudType>::g() const
{
    return g_;
}


template<class CloudType>
inline Foam::scalar Foam::FrontTrackingCloud<CloudType>::pAmbient() const
{
    return pAmbient_;
}


template<class CloudType>
inline Foam::scalar& Foam::FrontTrackingCloud<CloudType>::pAmbient()
{
    return pAmbient_;
}

/*
template<class CloudType>
//inline const typename CloudType::parcelType::forceType&
inline const typename Foam::FrontTrackingCloud<CloudType>::forceType&
Foam::FrontTrackingCloud<CloudType>::forces() const
{
    return forces_;
}


template<class CloudType>
inline typename Foam::FrontTrackingCloud<CloudType>::forceType&
Foam::FrontTrackingCloud<CloudType>::forces()
{
    return forces_;
}
*/

/*
template<class CloudType>
inline typename Foam::FrontTrackingCloud<CloudType>::functionType&
Foam::FrontTrackingCloud<CloudType>::functions()
{
    return functions_;
}
*/

/*
template<class CloudType>
inline const Foam::InjectionModelList<Foam::FrontTrackingCloud<CloudType> >&
Foam::FrontTrackingCloud<CloudType>::injectors() const
{
    return injectors_;
}


template<class CloudType>
inline Foam::InjectionModelList<Foam::FrontTrackingCloud<CloudType> >&
Foam::FrontTrackingCloud<CloudType>::injectors()
{
    return injectors_;
}


template<class CloudType>
inline const Foam::DispersionModel<Foam::FrontTrackingCloud<CloudType> >&
Foam::FrontTrackingCloud<CloudType>::dispersion() const
{
    return dispersionModel_;
}


template<class CloudType>
inline Foam::DispersionModel<Foam::FrontTrackingCloud<CloudType> >&
Foam::FrontTrackingCloud<CloudType>::dispersion()
{
    return dispersionModel_();
}
*/

template<class CloudType>
inline const Foam::FrontAnalyzingModel<Foam::FrontTrackingCloud<CloudType> >&
Foam::FrontTrackingCloud<CloudType>::FrontAnalyzing() const
{
    return FrontAnalyzingModel_;
}


template<class CloudType>
inline Foam::FrontAnalyzingModel<Foam::FrontTrackingCloud<CloudType> >&
Foam::FrontTrackingCloud<CloudType>::FrontAnalyzing()
{
    return FrontAnalyzingModel_();
}


/*
template<class CloudType>
inline const Foam::StochasticCollisionModel<Foam::FrontTrackingCloud<CloudType> >&
Foam::FrontTrackingCloud<CloudType>::stochasticCollision() const
{
    return stochasticCollisionModel_();
}


template<class CloudType>
inline Foam::StochasticCollisionModel<Foam::FrontTrackingCloud<CloudType> >&
Foam::FrontTrackingCloud<CloudType>::stochasticCollision()
{
    return stochasticCollisionModel_();
}


template<class CloudType>
inline const Foam::SurfaceFilmModel<Foam::FrontTrackingCloud<CloudType> >&
Foam::FrontTrackingCloud<CloudType>::surfaceFilm() const
{
    return surfaceFilmModel_();
}


template<class CloudType>
inline Foam::SurfaceFilmModel<Foam::FrontTrackingCloud<CloudType> >&
Foam::FrontTrackingCloud<CloudType>::surfaceFilm()
{
    return surfaceFilmModel_();
}


template<class CloudType> //added ttw
inline const Foam::ModulationModel<Foam::FrontTrackingCloud<CloudType> >&
Foam::FrontTrackingCloud<CloudType>::modulation() const
{
    return modulationModel_;
}


template<class CloudType>
inline Foam::ModulationModel<Foam::FrontTrackingCloud<CloudType> >&
Foam::FrontTrackingCloud<CloudType>::modulation()
{
    return modulationModel_();
}
*/

/*

template<class CloudType> //added1
inline const bool&
Foam::FrontTrackingCloud<CloudType>::preMotStat1Comp() const
{
    return preMotStat1Comp_;
}


template<class CloudType>
inline bool&
Foam::FrontTrackingCloud<CloudType>::preMotStat1Comp()
{
    return preMotStat1Comp_;
}


template<class CloudType>
inline const bool&
Foam::FrontTrackingCloud<CloudType>::preColStat1Comp() const
{
    return preColStat1Comp_;
}


template<class CloudType>
inline bool&
Foam::FrontTrackingCloud<CloudType>::preColStat1Comp()
{
    return preColStat1Comp_;
}


template<class CloudType> //added1
inline const bool&
Foam::FrontTrackingCloud<CloudType>::preMotStat2Comp() const
{
    return preMotStat2Comp_;
}


template<class CloudType>
inline bool&
Foam::FrontTrackingCloud<CloudType>::preMotStat2Comp()
{
    return preMotStat2Comp_;
}


template<class CloudType>
inline const bool&
Foam::FrontTrackingCloud<CloudType>::preColStat2Comp() const
{
    return preColStat2Comp_;
}


template<class CloudType>
inline bool&
Foam::FrontTrackingCloud<CloudType>::preColStat2Comp()
{
    return preColStat2Comp_;
}
*/

/*

template<class CloudType>
inline const Foam::vectorIntegrationScheme&
Foam::FrontTrackingCloud<CloudType>::UIntegrator() const
{
    return UIntegrator_;
}
*/

template<class CloudType>
inline Foam::label Foam::FrontTrackingCloud<CloudType>::nParcels() const
{
    return this->size();
}


template<class CloudType>
inline Foam::scalar Foam::FrontTrackingCloud<CloudType>::massInSystem() const
{
    scalar sysMass = 0.0;
    forAllConstIter(typename FrontTrackingCloud<CloudType>, *this, iter)
    {
         const parcelType& p = iter();
         sysMass += p.nParticle()*p.mass();
    }

    return sysMass;
}

/*
template<class CloudType>
inline Foam::vector
Foam::FrontTrackingCloud<CloudType>::linearMomentumOfSystem() const
{
    vector linearMomentum(vector::zero);

    forAllConstIter(typename FrontTrackingCloud<CloudType>, *this, iter)
    {
        const parcelType& p = iter();

        linearMomentum += p.nParticle()*p.mass()*p.U();
    }

    return linearMomentum;
}


template<class CloudType>
inline Foam::scalar
Foam::FrontTrackingCloud<CloudType>::linearKineticEnergyOfSystem() const
{
    scalar linearKineticEnergy = 0.0;

    forAllConstIter(typename FrontTrackingCloud<CloudType>, *this, iter)
    {
        const parcelType& p = iter();

        linearKineticEnergy += p.nParticle()*0.5*p.mass()*(p.U() & p.U());
    }

    return linearKineticEnergy;
}


template<class CloudType>
inline Foam::scalar Foam::FrontTrackingCloud<CloudType>::Dij
(
    const label i,
    const label j
) const
{
    scalar si = 0.0;
    scalar sj = 0.0;
    forAllConstIter(typename FrontTrackingCloud<CloudType>, *this, iter)
    {
        const parcelType& p = iter();
        si += p.nParticle()*pow(p.d(), i);
        sj += p.nParticle()*pow(p.d(), j);
    }

    reduce(si, sumOp<scalar>());
    reduce(sj, sumOp<scalar>());
    sj = max(sj, VSMALL);

    return si/sj;
}


template<class CloudType>
inline Foam::scalar Foam::FrontTrackingCloud<CloudType>::Dmax() const
{
    scalar d = -GREAT;
    forAllConstIter(typename FrontTrackingCloud<CloudType>, *this, iter)
    {
        const parcelType& p = iter();
        d = max(d, p.d());
    }

    reduce(d, maxOp<scalar>());

    return max(0.0, d);
}

*/

/*
template<class CloudType>
inline Foam::scalar Foam::FrontTrackingCloud<CloudType>::penetration
(
    const scalar fraction
) const
{
    if ((fraction < 0) || (fraction > 1))
    {
        FatalErrorIn
        (
            "inline Foam::scalar Foam::FrontTrackingCloud<CloudType>::penetration"
            "("
                "const scalar"
            ") const"
        )
            << "fraction should be in the range 0 < fraction < 1"
            << exit(FatalError);
    }

    scalar distance = 0.0;

    const label nParcel = this->size();
    globalIndex globalParcels(nParcel);
    const label nParcelSum = globalParcels.size();

    if (nParcelSum == 0)
    {
        return distance;
    }

    // lists of parcels mass and distance from initial injection point
    List<List<scalar> > procMass(Pstream::nProcs());
    List<List<scalar> > procDist(Pstream::nProcs());

    List<scalar>& mass = procMass[Pstream::myProcNo()];
    List<scalar>& dist = procDist[Pstream::myProcNo()];

    mass.setSize(nParcel);
    dist.setSize(nParcel);

    label i = 0;
    scalar mSum = 0.0;
    forAllConstIter(typename FrontTrackingCloud<CloudType>, *this, iter)
    {
        const parcelType& p = iter();
        scalar m = p.nParticle()*p.mass();
        scalar d = mag(p.position() - p.position0());
        mSum += m;

        mass[i] = m;
        dist[i] = d;

        i++;
    }

    // calculate total mass across all processors
    reduce(mSum, sumOp<scalar>());
    Pstream::gatherList(procMass);
    Pstream::gatherList(procDist);

    if (Pstream::master())
    {
        // flatten the mass lists
        List<scalar> allMass(nParcelSum, 0.0);
        SortableList<scalar> allDist(nParcelSum, 0.0);
        for (label procI = 0; procI < Pstream::nProcs(); procI++)
        {
            SubList<scalar>
            (
                allMass,
                globalParcels.localSize(procI),
                globalParcels.offset(procI)
            ).assign(procMass[procI]);

            // flatten the distance list
            SubList<scalar>
            (
                allDist,
                globalParcels.localSize(procI),
                globalParcels.offset(procI)
            ).assign(procDist[procI]);
        }

        // sort allDist distances into ascending order
        // note: allMass masses are left unsorted
        allDist.sort();

        if (nParcelSum > 1)
        {
            const scalar mLimit = fraction*mSum;
            const labelList& indices = allDist.indices();

            if (mLimit > (mSum - allMass[indices.last()]))
            {
                distance = allDist.last();
            }
            else
            {
                // assuming that 'fraction' is generally closer to 1 than 0,
                // loop through in reverse distance order
                const scalar mThreshold = (1.0 - fraction)*mSum;
                scalar mCurrent = 0.0;
                label i0 = 0;

                forAllReverse(indices, i)
                {
                    label indI = indices[i];

                    mCurrent += allMass[indI];

                    if (mCurrent > mThreshold)
                    {
                        i0 = i;
                        break;
                    }
                }

                if (i0 == indices.size() - 1)
                {
                    distance = allDist.last();
                }
                else
                {
                    // linearly interpolate to determine distance
                    scalar alpha = (mCurrent - mThreshold)/allMass[indices[i0]];
                    distance =
                        allDist[i0] + alpha*(allDist[i0+1] - allDist[i0]);
                }
            }
        }
        else
        {
            distance = allDist.first();
        }
    }

    Pstream::scatter(distance);

    return distance;
}
*/

template<class CloudType>
inline Foam::cachedRandom& Foam::FrontTrackingCloud<CloudType>::rndGen()
{
    return rndGen_;
}


template<class CloudType>
inline Foam::List<Foam::DynamicList<typename CloudType::particleType*> >&
Foam::FrontTrackingCloud<CloudType>::cellOccupancy()
{
    if (cellOccupancyPtr_.empty())
    {
        buildCellOccupancy();
    }

    return cellOccupancyPtr_();
}


template<class CloudType>
inline const Foam::scalarField&
Foam::FrontTrackingCloud<CloudType>::cellLengthScale() const
{
    return cellLengthScale_;
}

/*
template<class CloudType>
inline Foam::DimensionedField<Foam::vector, Foam::volMesh>&
Foam::FrontTrackingCloud<CloudType>::UTrans()
{
    return UTrans_();
}


template<class CloudType>
inline const Foam::DimensionedField<Foam::vector, Foam::volMesh>&
Foam::FrontTrackingCloud<CloudType>::UTrans() const
{
    return UTrans_();
}


template<class CloudType>
inline Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::FrontTrackingCloud<CloudType>::UCoeff()
{
    return UCoeff_();
}


template<class CloudType>
inline const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::FrontTrackingCloud<CloudType>::UCoeff() const
{
    return UCoeff_();
}


template<class CloudType>
inline Foam::DimensionedField<Foam::symmTensor, Foam::volMesh>&
Foam::FrontTrackingCloud<CloudType>::RTrans()
{
    return RTrans_();
}


template<class CloudType>
inline const Foam::DimensionedField<Foam::symmTensor, Foam::volMesh>&
Foam::FrontTrackingCloud<CloudType>::RTrans() const
{
    return RTrans_();
}
*/

//Jafari added.
template<class CloudType>
inline Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::FrontTrackingCloud<CloudType>::densityFromFT()
{
    return densityFromFT_();
}

//Jafari added.
template<class CloudType>
inline const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::FrontTrackingCloud<CloudType>::densityFromFT() const
{
    return densityFromFT_();
}


//Jafari added.
template<class CloudType>
inline Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::FrontTrackingCloud<CloudType>::viscosityFromFT()
{
    return viscosityFromFT_();
}

//Jafari added.
template<class CloudType>
inline const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::FrontTrackingCloud<CloudType>::viscosityFromFT() const
{
    return viscosityFromFT_();
}


//Jafari added.
template<class CloudType>
inline Foam::DimensionedField<Foam::vector, Foam::volMesh>&
Foam::FrontTrackingCloud<CloudType>::sTensionForceFromFT()
{
    return sTensionForceFromFT_();
}

//Jafari added.
template<class CloudType>
inline const Foam::DimensionedField<Foam::vector, Foam::volMesh>&
Foam::FrontTrackingCloud<CloudType>::sTensionForceFromFT() const
{
    return sTensionForceFromFT_();
}

//Jafari added.
template<class CloudType>
inline Foam::DimensionedField<Foam::vector, Foam::volMesh>&
Foam::FrontTrackingCloud<CloudType>::dPressureFromFT()
{
    return dPressureFromFT_();
}

//Jafari added.
template<class CloudType>
inline const Foam::DimensionedField<Foam::vector, Foam::volMesh>&
Foam::FrontTrackingCloud<CloudType>::dPressureFromFT() const
{
    return dPressureFromFT_();
}

//Jafari added.
template<class CloudType>
inline Foam::DimensionedField<Foam::vector, Foam::volMesh>&
Foam::FrontTrackingCloud<CloudType>::pressureJumpAtTheInterfaceFromFT()
{
    return pressureJumpAtTheInterfaceFromFT_();
}

//Jafari added.
template<class CloudType>
inline const Foam::DimensionedField<Foam::vector, Foam::volMesh>&
Foam::FrontTrackingCloud<CloudType>::pressureJumpAtTheInterfaceFromFT() const
{
    return pressureJumpAtTheInterfaceFromFT_();
}

/*
template<class CloudType>
inline Foam::tmp<Foam::fvVectorMatrix>
Foam::FrontTrackingCloud<CloudType>::SU(volVectorField& U) const
{
    if (debug)
    {
        Info<< "UTrans min/max = " << min(UTrans()).value() << ", "
            << max(UTrans()).value() << nl
            << "UCoeff min/max = " << min(UCoeff()).value() << ", "
            << max(UCoeff()).value() << endl;
    }

    if (solution_.coupled())
    {
        if (solution_.semiImplicit("U"))
        {
            const DimensionedField<scalar, volMesh>
                Vdt(mesh_.V()*this->db().time().deltaT());

            return UTrans()/Vdt - fvm::Sp(UCoeff()/Vdt, U) + UCoeff()/Vdt*U;
        }
        else
        {
            tmp<fvVectorMatrix> tfvm(new fvVectorMatrix(U, dimForce));
            fvVectorMatrix& fvm = tfvm();

            fvm.source() = -UTrans()/(this->db().time().deltaT());

            return tfvm;
        }
    }

    return tmp<fvVectorMatrix>(new fvVectorMatrix(U, dimForce));
}
*/


/*template<class CloudType> //added ttw
inline Foam::tmp<Foam::fvSymmTensorMatrix>
Foam::FrontTrackingCloud<CloudType>::SR(volSymmTensorField& R) const
{
    if (debug)
    {
        Info<< "RTrans min/max = " << min(RTrans()).value() << ", "
            << max(RTrans()).value() << endl;
    }

	if (solution_.turbulenceCoupling())
    {
        if (solution_.semiImplicit("U"))
        {
            const DimensionedField<scalar, volMesh>
                Vdt(mesh_.V()*this->db().time().deltaT());

            return RTrans()/Vdt;
        }
        else
        {
            tmp<fvSymmTensorMatrix> tfvm(new fvSymmTensorMatrix(R, dimForce*dimVelocity));
		    fvSymmTensorMatrix& fvm = tfvm();

		    fvm.source() = -RTrans()/(this->db().time().deltaT());

		    return tfvm;
        }
    }

    return tmp<fvSymmTensorMatrix>(new fvSymmTensorMatrix(R, dimForce*dimVelocity));
}*/

/*
template<class CloudType> //added ttw
inline Foam::tmp<Foam::volSymmTensorField>
Foam::FrontTrackingCloud<CloudType>::SR() const
{
    if (debug)
    {
        Info<< "RTrans min/max = " << min(RTrans()).value() << ", "
            << max(RTrans()).value() << endl;
    }

	tmp<volSymmTensorField> tSR
    (
        new volSymmTensorField
        (
            IOobject
            (
                this->name() + ":SR",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedSymmTensor("zero", dimForce*dimVelocity/dimVolume, symmTensor::zero),
            zeroGradientFvPatchSymmTensorField::typeName
        )
    );
    volSymmTensorField& SR = tSR();

	if (solution_.turbulenceCoupling()) SR.internalField() = RTrans()/(mesh_.V()*this->db().time().deltaT());

    return tSR;
}
*/

//Jafari added.
template<class CloudType>
inline const Foam::tmp<Foam::volScalarField>
Foam::FrontTrackingCloud<CloudType>::densityIndicator() const
{
    tmp<volScalarField> tdensityIndicator
    (
        new volScalarField
        (
            IOobject
            (
                this->name() + ":densityIndicator",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimDensity, 0.0)//,
            //zeroGradientFvPatchScalarField::typeName
        )
    );

    volScalarField& densityIndicator = tdensityIndicator();

    densityIndicator.internalField() = densityFromFT();
    densityIndicator.correctBoundaryConditions();

    return tdensityIndicator;
}

//Jafari added.
template<class CloudType>
inline const Foam::tmp<Foam::volScalarField>
Foam::FrontTrackingCloud<CloudType>::viscosityIndicator() const
{
    tmp<volScalarField> tviscosityIndicator
    (
        new volScalarField
        (
            IOobject
            (
                this->name() + ":viscosityIndicator",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimMass/dimLength/dimTime, 0.0)//,
            //zeroGradientFvPatchScalarField::typeName
        )
    );

    volScalarField& viscosityIndicator = tviscosityIndicator();

    viscosityIndicator.internalField() = viscosityFromFT();
    viscosityIndicator.correctBoundaryConditions();

    return tviscosityIndicator;
}

//Jafari added.
template<class CloudType>
inline const Foam::tmp<Foam::volVectorField>
Foam::FrontTrackingCloud<CloudType>::sTensionForce() const
{
    tmp<volVectorField> tsTensionForce
    (
        new volVectorField
        (
            IOobject
            (
                this->name() + ":sTensionForce",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedVector("vector::zero", dimForce/dimVolume,  vector::zero)
        )
    );

    volVectorField& sTensionForce = tsTensionForce();

    sTensionForce.internalField() = sTensionForceFromFT();///mesh_.V();
    sTensionForce.correctBoundaryConditions();

    return tsTensionForce;
}

//Jafari added.
template<class CloudType>
inline const Foam::tmp<Foam::volVectorField>
Foam::FrontTrackingCloud<CloudType>::dPressure() const
{
    tmp<volVectorField> tdPressure
    (
        new volVectorField
        (
            IOobject
            (
                this->name() + ":dPressure",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedVector("vector::zero", dimForce/dimVolume,  vector::zero)
        )
    );

    volVectorField& dPressure = tdPressure();

    dPressure.internalField() = dPressureFromFT();///mesh_.V();
    dPressure.correctBoundaryConditions();

    return tdPressure;
}

//Jafari added.
template<class CloudType>
inline const Foam::tmp<Foam::volVectorField>
Foam::FrontTrackingCloud<CloudType>::pressureJumpAtTheInterface() const
{
    tmp<volVectorField> tpressureJumpAtTheInterface
    (
        new volVectorField
        (
            IOobject
            (
                this->name() + ":pressureJumpAtTheInterface",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedVector("vector::zero", dimForce/dimVolume,  vector::zero)
        )
    );

    volVectorField& pressureJumpAtTheInterface = tpressureJumpAtTheInterface();

    pressureJumpAtTheInterface.internalField() = pressureJumpAtTheInterfaceFromFT();///mesh_.V();
    pressureJumpAtTheInterface.correctBoundaryConditions();

    return tpressureJumpAtTheInterface;
}

template<class CloudType>
inline const Foam::tmp<Foam::volScalarField>
Foam::FrontTrackingCloud<CloudType>::vDotSweep() const
{
    tmp<volScalarField> tvDotSweep
    (
        new volScalarField
        (
            IOobject
            (
                this->name() + ":vDotSweep",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimless/dimTime, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    volScalarField& vDotSweep = tvDotSweep();
    forAllConstIter(typename FrontTrackingCloud<CloudType>, *this, iter)
    {
        const parcelType& p = iter();
        const label cellI = p.cell();

        vDotSweep[cellI] += p.nParticle()*p.areaP()*mag(p.U() - U_[cellI]);
    }

    vDotSweep.internalField() /= mesh_.V();
    vDotSweep.correctBoundaryConditions();

    return tvDotSweep;
}


template<class CloudType>
inline const Foam::tmp<Foam::volScalarField>
Foam::FrontTrackingCloud<CloudType>::theta() const
{
    tmp<volScalarField> ttheta
    (
        new volScalarField
        (
            IOobject
            (
                this->name() + ":theta",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimless, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    volScalarField& theta = ttheta();
    forAllConstIter(typename FrontTrackingCloud<CloudType>, *this, iter)
    {
        const parcelType& p = iter();
        const label cellI = p.cell();

        theta[cellI] += p.nParticle()*p.volume();
    }

    theta.internalField() /= mesh_.V();
    theta.correctBoundaryConditions();

    return ttheta;
}


template<class CloudType>
inline const Foam::tmp<Foam::volScalarField>
Foam::FrontTrackingCloud<CloudType>::alpha() const
{
    tmp<volScalarField> talpha
    (
        new volScalarField
        (
            IOobject
            (
                this->name() + ":alpha",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimless, 0.0)
        )
    );

    scalarField& alpha = talpha().internalField();
    forAllConstIter(typename FrontTrackingCloud<CloudType>, *this, iter)
    {
        const parcelType& p = iter();
        const label cellI = p.cell();

        alpha[cellI] += p.nParticle()*p.mass();
    }

    alpha /= (mesh_.V()*rho_);

    return talpha;
}


template<class CloudType>
inline const Foam::tmp<Foam::volScalarField>
Foam::FrontTrackingCloud<CloudType>::rhoEff() const
{
    tmp<volScalarField> trhoEff
    (
        new volScalarField
        (
            IOobject
            (
                this->name() + ":rhoEff",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimDensity, 0.0)
        )
    );

    scalarField& rhoEff = trhoEff().internalField();
    forAllConstIter(typename FrontTrackingCloud<CloudType>, *this, iter)
    {
        const parcelType& p = iter();
        const label cellI = p.cell();

        rhoEff[cellI] += p.nParticle()*p.mass();
    }

    rhoEff /= mesh_.V();

    return trhoEff;
}


// ************************************************************************* //
