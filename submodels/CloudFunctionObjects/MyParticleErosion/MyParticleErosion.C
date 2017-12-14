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

#include "MyParticleErosion.H"

// * * * * * * * * * * * * * Protectd Member Functions * * * * * * * * * * * //

template<class CloudType>
Foam::label Foam::MyParticleErosion<CloudType>::applyToPatch
(
    const label globalPatchI
) const
{
    forAll(patchIDs_, i)
    {
        if (patchIDs_[i] == globalPatchI)
        {
            return i;
        }
    }

    return -1;
}


template<class CloudType>
void Foam::MyParticleErosion<CloudType>::write()
{
    if (QPtr_.valid())
    {
        QPtr_->write();
    }
    else
    {
        FatalErrorIn("void Foam::MyParticleErosion<CloudType>::write()")
            << "QPtr not valid" << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::MyParticleErosion<CloudType>::MyParticleErosion
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName
)
:
    CloudFunctionObject<CloudType>(dict, owner, modelName, typeName),
    QPtr_(NULL),
    patchIDs_(),
    ner_(readScalar(this->coeffDict().lookup("ner"))),
    material_(readScalar(this->coeffDict().lookup("material")))
{
    const wordList allPatchNames = owner.mesh().boundaryMesh().names();
    wordList patchName(this->coeffDict().lookup("patches"));

    labelHashSet uniquePatchIDs;
    forAllReverse(patchName, i)
    {
        labelList patchIDs = findStrings(patchName[i], allPatchNames);

        if (patchIDs.empty())
        {
            WarningIn
            (
                "Foam::MyParticleErosion<CloudType>::MyParticleErosion"
                "("
                    "const dictionary&, "
                    "CloudType& "
                ")"
            )   << "Cannot find any patch names matching " << patchName[i]
                << endl;
        }

        uniquePatchIDs.insert(patchIDs);
    }

    patchIDs_ = uniquePatchIDs.toc();

    // trigger ther creation of the Q field
    preEvolve();
}


template<class CloudType>
Foam::MyParticleErosion<CloudType>::MyParticleErosion
(
    const MyParticleErosion<CloudType>& pe
)
:
    CloudFunctionObject<CloudType>(pe),
    QPtr_(NULL),
    patchIDs_(pe.patchIDs_),
    ner_(pe.ner_),
    material_(pe.material_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::MyParticleErosion<CloudType>::~MyParticleErosion()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::MyParticleErosion<CloudType>::preEvolve()
{
    if (QPtr_.valid())
    {
        QPtr_->internalField() = 0.0;
    }
    else
    {
        const fvMesh& mesh = this->owner().mesh();

        QPtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    this->owner().name() + "Q",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar("zero", dimVolume, 0.0)
            )
        );
    }
}


template<class CloudType>
void Foam::MyParticleErosion<CloudType>::postPatch
(
    const parcelType& p,
    const polyPatch& pp,
    const scalar trackFraction,
    const tetIndices& tetIs,
    bool&
)
{
    const label patchI = pp.index();

    const label localPatchI = applyToPatch(patchI);

    if (localPatchI != -1)
    {
        vector nw;
        vector Up;

        // patch-normal direction
        this->owner().patchData(p, pp, trackFraction, tetIs, nw, Up);

        // particle velocity reletive to patch
        const vector& U = p.U() - Up;

        // quick reject if particle travelling away from the patch
        if ((nw & U) < 0)
        {
            return;
        }

        const scalar magU = mag(U);
        const vector Udir = U/magU;

        // determine impact angle, alpha
        const scalar alpha = mathematical::pi/2.0 - acos(nw & Udir);

	const scalar ER1=material_/pow(p.d(),(1.-ner_)/4.);  
	const scalar ER2=pow(p.mass(),(1.+3.*(1-ner_)/4.0));
	const scalar ER3=pow(magU,(2.+3.*(1-ner_)/2.0));
	const scalar ER4=pow(cos(alpha),2.0);
	const scalar ER5=pow(sin(alpha),(3.*(1.-ner_)/2.0));
	const scalar erosion=365.*24*3600.*ER1*ER2*ER3*ER4*ER5*p.nParticle();

        const label patchFaceI = pp.whichFace(p.face());
        scalar& Q = QPtr_->boundaryField()[patchI][patchFaceI];

	Q += erosion;

        /*
        if (tan(alpha) < 1/6.0)
        {
            Q += coeff*(sin(2.0*alpha) - 6.0/1*sqr(sin(alpha)));
        }
        else
        {
            Q += coeff*(1*sqr(cos(alpha))/6.0);
        }
        */

    }
}


// ************************************************************************* //
