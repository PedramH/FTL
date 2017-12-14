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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
inline const Foam::KinematicRotCloud<CloudType>&
Foam::KinematicRotCloud<CloudType>::cloudCopy() const
{
    return cloudCopyPtr_();
}


template<class CloudType>
inline const typename CloudType::particleType::constantProperties&
Foam::KinematicRotCloud<CloudType>::constProps() const
{
    return constProps_;
}


template<class CloudType>
inline typename CloudType::particleType::constantProperties&
Foam::KinematicRotCloud<CloudType>::constProps()
{
    return constProps_;
}


template<class CloudType>
inline const typename Foam::KinematicRotCloud<CloudType>::torqueType&
Foam::KinematicRotCloud<CloudType>::torques() const
{
    return torques_;
}


template<class CloudType>
inline typename Foam::KinematicRotCloud<CloudType>::torqueType&
Foam::KinematicRotCloud<CloudType>::torques()
{
    return torques_;
}


/*template<class CloudType>
inline const Foam::DeterminesticCollisionModel<Foam::KinematicRotCloud<CloudType> >&
Foam::KinematicRotCloud<CloudType>::determinesticCollision() const
{
    return determinesticCollisionModel_();
}


template<class CloudType>
inline Foam::DeterminesticCollisionModel<Foam::KinematicRotCloud<CloudType> >&
Foam::KinematicRotCloud<CloudType>::determinesticCollision()
{
    return determinesticCollisionModel_();
}*/

//OmegaC memory
/*template<class CloudType>
inline const Foam::volVectorField& Foam::KinematicRotCloud<CloudType>::curlU() const
{
	if (!curlUPtr_.valid())
    {
        FatalErrorIn
        (
            "inline const Foam::volVectorField&"
            "Foam::KinematicRotCloud<CloudType>::curlU() const"
        )   << "Carrier phase curlUc Field not set"
            << abort(FatalError);
    }

    return curlUPtr_();
}*/


/*template<class CloudType>
inline const bool&
Foam::KinematicRotCloud<CloudType>::stat2Comp() const
{
    return stat2Comp_;
}


template<class CloudType>
inline bool&
Foam::KinematicRotCloud<CloudType>::stat2Comp()
{
    return stat2Comp_;
}*/


template<class CloudType>
inline const Foam::vectorIntegrationScheme&
Foam::KinematicRotCloud<CloudType>::OmegaIntegrator() const
{
    return OmegaIntegrator_;
}


template<class CloudType>
inline Foam::scalar Foam::KinematicRotCloud<CloudType>::Omegamax() const
{
    scalar T = -GREAT;
    scalar n = 0;
    forAllConstIter(typename KinematicRotCloud<CloudType>, *this, iter)
    {
        const parcelType& p = iter();
        T = max(T, mag(p.Omega()));
        n++;
    }

    reduce(T, maxOp<scalar>());
    reduce(n, sumOp<label>());

    if (n > 0)
    {
        return T;
    }
    else
    {
        return 0.0;
    }
}


template<class CloudType>
inline Foam::scalar Foam::KinematicRotCloud<CloudType>::Omegamin() const
{
    scalar T = GREAT;
    scalar n = 0;
    forAllConstIter(typename KinematicRotCloud<CloudType>, *this, iter)
    {
        const parcelType& p = iter();
        T = min(T, mag(p.Omega()));
        n++;
    }

    reduce(T, minOp<scalar>());
    reduce(n, sumOp<label>());

    if (n > 0)
    {
        return T;
    }
    else
    {
        return 0.0;
    }
}


template<class CloudType>
inline Foam::DimensionedField<Foam::vector, Foam::volMesh>&
Foam::KinematicRotCloud<CloudType>::OmegaTrans()
{
    return OmegaTrans_();
}


template<class CloudType>
inline const Foam::DimensionedField<Foam::vector, Foam::volMesh>&
Foam::KinematicRotCloud<CloudType>::OmegaTrans() const
{
    return OmegaTrans_();
}

/*
template<class CloudType>
inline Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::KinematicRotCloud<CloudType>::OmegaCoeff()
{
    return OmegaCoeff_();
}


template<class CloudType>
inline const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::KinematicRotCloud<CloudType>::OmegaCoeff() const
{
    return OmegaCoeff_();
}
*/


/*template<class CloudType>
inline const Foam::tmp<Foam::volScalarField>
Foam::KinematicRotCloud<CloudType>::NP() const
{
    tmp<volScalarField> tNP
    (
        new volScalarField
        (
            IOobject
            (
                this->name() + ":NP",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", 1.0/dimVolume, 0.0)
        )
    );

    scalarField& NP = tNP().internalField();
    forAllConstIter(typename KinematicCloud<CloudType>, *this, iter)
    {
        const parcelType& p = iter();
        const label cellI = p.cell();

        NP[cellI] += p.nParticle();
    }

    //NP /= mesh_.V();

    return tNP;
}


template<class CloudType>
inline const Foam::tmp<Foam::volScalarField>
Foam::KinematicRotCloud<CloudType>::MP() const
{
    tmp<volScalarField> tMP
    (
        new volScalarField
        (
            IOobject
            (
                this->name() + ":MP",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", 1.0/dimVolume, 0.0)
        )
    );

    scalarField& MP = tMP().internalField();
    forAllConstIter(typename KinematicCloud<CloudType>, *this, iter)
    {
        const parcelType& p = iter();
        const label cellI = p.cell();

        MP[cellI] += p.nParticle()*p.mass();
    }

    //MP /= mesh_.V();

    return tMP;
}


template<class CloudType>
inline const Foam::tmp<Foam::volVectorField>
Foam::KinematicRotCloud<CloudType>::UP() const
{
    tmp<volVectorField> tUP
    (
        new volVectorField
        (
            IOobject
            (
                this->name() + ":UP",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedVelocity("zero", dimVelocity, vector::zero)
        )
    );

    vectorField& UP = tUP().internalField();
    forAllConstIter(typename KinematicCloud<CloudType>, *this, iter)
    {
        const parcelType& p = iter();
        const label cellI = p.cell();

		const scalar m = p.nParticle()*p.mass();

        UP[cellI] += m*p.U();
    }

    UP /= max(MP(), SMALL); //MP();

    return tUP;
}


template<class CloudType>
inline const Foam::tmp<Foam::volVectorField>
Foam::KinematicRotCloud<CloudType>::UPUP() const
{
    tmp<volVectorField> tUPUP
    (
        new volVectorField
        (
            IOobject
            (
                this->name() + ":UPUP",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedVelocity("zero", dimVelocity*dimVelocity, vector::zero)
        )
    );

    vectorField& UPUP = tUPUP().internalField();
    forAllConstIter(typename KinematicCloud<CloudType>, *this, iter)
    {
        const parcelType& p = iter();
        const label cellI = p.cell();

		const scalar m = p.nParticle()*p.mass();

        UPUP[cellI] += m*p.U()*p.U();
    }

    UPUP /= max(MP(), SMALL); //MP();

    return tUPUP;
}


template<class CloudType>
inline const Foam::tmp<Foam::volScalarField>
Foam::KinematicRotCloud<CloudType>::d10() const
{
    tmp<volScalarField> td10
    (
        new volScalarField
        (
            IOobject
            (
                this->name() + ":d10",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimLength, 0.0)
        )
    );

    scalarField& d10 = td10().internalField();
    forAllConstIter(typename KinematicCloud<CloudType>, *this, iter)
    {
        const parcelType& p = iter();
        const label cellI = p.cell();

        d10[cellI] += p.nParticle()*p.d();
    }

    d10 /= max(NP(), SMALL); //NP();

    return td10;
}
*/

// ************************************************************************* //
