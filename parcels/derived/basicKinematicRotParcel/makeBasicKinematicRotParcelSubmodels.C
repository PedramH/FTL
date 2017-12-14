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

#include "basicKinematicRotCloud.H"

#include "makeParcelCloudFunctionObjects.H"

// Kinematic
#include "makeKinematicRotParcelForces.H" // KinematicRot variant
#include "makeKinematicRotParcelTorques.H" 
#include "makeParcelDispersionModels.H"
#include "makeParcelInjectionModels.H"
#include "makeParcelPatchInteractionModels.H"
#include "makeKinematicRotParcelStochasticCollisionModels.H" // KinematicRot variant
#include "makeParcelSurfaceFilmModels.H"
#include "makeParcelModulationModels.H" //added ttw

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makeParcelCloudFunctionObjects(basicKinematicRotCloud);

    // Kinematic sub-models
    makeKinematicRotParcelForces(basicKinematicRotCloud); // KinematicRot variant
	makeKinematicRotParcelTorques(basicKinematicRotCloud);
    makeParcelDispersionModels(basicKinematicRotCloud);
    makeParcelInjectionModels(basicKinematicRotCloud);
    makeParcelPatchInteractionModels(basicKinematicRotCloud);
    makeKinematicRotParcelStochasticCollisionModels(basicKinematicRotCloud); // KinematicRot variant
    makeParcelSurfaceFilmModels(basicKinematicRotCloud);
	makeParcelModulationModels(basicKinematicRotCloud); //added ttw
}


// ************************************************************************* //
