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

#include "basicElectrokineticMPPICCloud.H"

#include "makeParcelCloudFunctionObjects.H"

// Kinematic
//#include "makeElectrokineticMPPICParcelForces.H" // ElectrokineticMPPIC variant 
#include "makeElectrokineticParcelForces.H" // Electrokinetic variant
#include "makeKinematicRotParcelTorques.H" 
#include "makeParcelDispersionModels.H"
#include "makeElectrokineticParcelInjectionModels.H" // Electrokinetic variant
#include "makeParcelPatchInteractionModels.H"
//#include "makeParcelStochasticCollisionModels.H"
#include "makeKinematicRotParcelStochasticCollisionModels.H" // KinematicRot variant
#include "makeParcelSurfaceFilmModels.H"
#include "makeParcelModulationModels.H" //added ttw

// MPPIC sub-models
#include "makeMPPICParcelDampingModels.H"
#include "makeMPPICParcelIsotropyModels.H"
#include "makeMPPICParcelPackingModels.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makeParcelCloudFunctionObjects(basicElectrokineticMPPICCloud);

    // Kinematic sub-models
	//makeElectrokineticMPPICParcelForces(basicElectrokineticMPPICCloud); // ElectrokineticMPPIC variant
    makeElectrokineticParcelForces(basicElectrokineticMPPICCloud); // Electrokinetic variant
	makeKinematicRotParcelTorques(basicElectrokineticMPPICCloud);
    makeParcelDispersionModels(basicElectrokineticMPPICCloud);
    makeElectrokineticParcelInjectionModels(basicElectrokineticMPPICCloud); // Electrokinetic variant
    makeParcelPatchInteractionModels(basicElectrokineticMPPICCloud);
    //makeParcelStochasticCollisionModels(basicElectrokineticMPPICCloud);  
	makeKinematicRotParcelStochasticCollisionModels(basicElectrokineticMPPICCloud); // KinematicRot variant
    makeParcelSurfaceFilmModels(basicElectrokineticMPPICCloud);
	makeParcelModulationModels(basicElectrokineticMPPICCloud); //added ttw

	// MPPIC sub-models
    makeMPPICParcelDampingModels(basicElectrokineticMPPICCloud);
    makeMPPICParcelIsotropyModels(basicElectrokineticMPPICCloud);
    makeMPPICParcelPackingModels(basicElectrokineticMPPICCloud);
}


// ************************************************************************* //
