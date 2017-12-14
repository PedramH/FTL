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

Class
    Foam::basicElectrokineticCollidingParcel

Description
    Definition of basic Electrokinetic colliding parcel

SourceFiles
    basicElectrokineticParcel.H 

\*---------------------------------------------------------------------------*/

#ifndef basicElectrokineticCollidingParcel_H
#define basicElectrokineticCollidingParcel_H

#include "contiguous.H"
#include "particle.H"
#include "KinematicParcel.H"
#include "KinematicRotParcel.H" 
#include "ElectrokineticParcel.H" 
#include "CollidingParcel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    typedef CollidingParcel<ElectrokineticParcel<KinematicRotParcel<KinematicParcel<particle> > > > //changed
        basicElectrokineticCollidingParcel;

    template<>
    inline bool contiguous<basicElectrokineticCollidingParcel>()
    {
        return true;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
