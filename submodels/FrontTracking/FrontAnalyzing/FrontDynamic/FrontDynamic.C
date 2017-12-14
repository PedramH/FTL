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

#include "FrontDynamic.H"
#include "meshTools.H"
#include "simpleMatrix.H"
using namespace Foam::constant;

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //
template<class CloudType>
void Foam::FrontDynamic<CloudType>::readInitialBubbles()
{
    const scalar time = this->owner().db().time().value();
    if ( time <= timeTobeSteady_ )
    {
        cloudCleaning();
        makingBubblesAtTheZeroTime();
        calculatingThePosInDomain();
        surfaceTensionForceDistribution();
        communicationsOfTheFrontAndEulerianGrid();
    }
    else if ( time > timeTobeSteady_ )
    {
        cloudCleaning();
        readDataOfTheFrontsAtTheCurrentTimeIntimeName();
        calculatingThePosInDomain();
        surfaceTensionForceDistribution();
        communicationsOfTheFrontAndEulerianGrid();
        mapFrontToParcel();
        this->owner().writeFields();
        //scalar startTime = this->owner().db().time().elapsedCpuTime();
        //scalar endTime = this->owner().db().time().elapsedCpuTime();
        //Info << "\n    The taken time for p insertion is " << (endTime - startTime) 
        //     << ' ' << "and the cloud size is " << ' ' << this->owner().size() << endl;
    }
}

template<class CloudType>
void Foam::FrontDynamic<CloudType>::dynamic()
{
   const scalar time = this->owner().db().time().value();
   const scalar deltaTime = this->owner().db().time().deltaTValue();

   printingAllImportantRunTimeData();

   if( time < timeTobeSteady_)
   {  
       Info << "\n----------> At the time: " << this->owner().db().time().value()
            << " the bubbles are fixed and have no motions to have steady state flow ....... " << endl;    
   }
   if( time == timeTobeSteady_ || ( time > timeTobeSteady_ && (time-deltaTime) <= timeTobeSteady_ ) )
   {  
       Info << "\n----------> At the time: " << this->owner().db().time().value()
            << " all points required are to be inserted ....... " << endl;
       mapFrontToParcel();
       this->owner().writeFields();
   }
   if( time > timeTobeSteady_)
   {  
       mapParcelToFront();
       cloudCleaning();
       processingThefrontMeshes();    
       calculatingThePosInDomain();
       surfaceTensionForceDistribution();
       communicationsOfTheFrontAndEulerianGrid();
       outputTimePostProcessing(); 
       mapFrontToParcel();
   }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
template<class CloudType>
Foam::FrontDynamic<CloudType>::FrontDynamic
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName
)
:
    FrontAnalyzingModel<CloudType>(dict, owner, modelName),
    ranGen_(label(0)),
    p_(16), 
    q_(16),
    r_(16),
    //mComponent_(label(4)),
    beta_(4),
    densityIndicator_
    (
        IOobject
        (
            "densityIndicator_",
            this->owner().db().time().timeName(),
            this->owner().mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->owner().mesh()
    ),
    viscosityIndicator_
    (
        IOobject
        (
            "viscosityIndicator_",
            this->owner().db().time().timeName(),
            this->owner().mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->owner().mesh()
    ),
    Indicator_
    (
        IOobject
        (
            "Indicator_",
            this->owner().db().time().timeName(),
            this->owner().mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->owner().mesh()
    ),
    cellsToRefine_
    (
        IOobject
        (
            "cellsToRefine_",
            this->owner().db().time().timeName(),
            this->owner().mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->owner().mesh(),
        dimensionedScalar("cellsToRefine", dimless, 0.0)
    )
{
    /*
    typename CloudType::parcelType::template
            TrackingData<CloudType> td_(this->owner()); 
    */

    timeOld_ = this->owner().db().time().value();//before main loop in solver.C, this is the start time ( equal to Old time).

    this->coeffDict().lookup("initialMakingFileOfTheBubbles") >> initialMakingFileOfTheBubbles_;
    this->coeffDict().lookup("periodicalOption") >> periodicalOption_;
    this->coeffDict().lookup("frontMeshInput") >> frontMeshInputOption_;
    this->coeffDict().lookup("undulationRemovalingOption") >> undulationRemovalingOption_;
    this->coeffDict().lookup("undulationRemovalingInterval") >> undulationRemovalingInterval_;
    this->coeffDict().lookup("undulationRemovalingType") >> undulationRemovalingType_;
    this->coeffDict().lookup("timeTobeSteady") >> timeTobeSteady_;
    this->coeffDict().lookup("smoothingAndSpreadingOption") >> smoothingAndSpreadingOption_;
    this->coeffDict().lookup("constructingTheIndicatorsOption") >> constructingTheIndicatorsOption_;
    this->coeffDict().lookup("surfaceTensionForceDistributionOption") >> surfaceTensionForceDistributionOption_;
    this->coeffDict().lookup("twoFluidFlow") >> twoFluidFlow_;
    this->coeffDict().lookup("searchInBox") >> searchInBox_;
    this->coeffDict().lookup("pressureJumpAtTheInterfaceOption") >> pressureJumpAtTheInterfaceOption_;

    // this commented part is important to use dictionaries
    /*
    const word& cloudName = this->owner().name();
    const fvMesh& mesh = this->owner().mesh();
    IOdictionary myProberDict
    		(
    		    IOobject
    		    (
    		        cloudName + "Properties",
    		        mesh.time().constant(),
    		        mesh,
    		        IOobject::MUST_READ_IF_MODIFIED,
    		        IOobject::NO_WRITE
    		    )
    		);
    */
    
    // Here, the dict is the same dictionary as subModels in (cloudName + "Properties") file.
    dictionary FrontDynamicCoeffsDict = dict.subOrEmptyDict("FrontDynamicCoeffs"); 
    dictionary periodicalOptionCoeffsDict = FrontDynamicCoeffsDict.subOrEmptyDict("periodicalOptionCoeffs");
    dictionary manuallyMeshInputCoeffsDict = FrontDynamicCoeffsDict.subOrEmptyDict("manuallyMeshInputCoeffs");
    dictionary finerCalculatedMeshInputCoeffsDict = FrontDynamicCoeffsDict.subOrEmptyDict("finerCalculatedMeshInputCoeffs");
    dictionary smoothingAndSpreadingManuallyCoeffs = FrontDynamicCoeffsDict.subOrEmptyDict("smoothingAndSpreadingManuallyCoeffs");      
    dictionary smoothingAndSpreadingCalculatedCoeffs = FrontDynamicCoeffsDict.subOrEmptyDict("smoothingAndSpreadingCalculatedCoeffs");
    dictionary constructingTheIndicatorsOptionCoeffs = FrontDynamicCoeffsDict.subOrEmptyDict("constructingTheIndicatorsOptionCoeffs");

    if( periodicalOption_ == "rotationalPeriodicalInCurvedDuct")
    {
        periodicalOptionCoeffsDict.lookup("thetaP") >> thetaP_; 
	//thetaP_ = readScalar(myDictL2.lookup("thetaP"));
        //myDictL2.readIfPresent("thetaP", thetaP_);
        periodicalOptionCoeffsDict.lookup("rCentre") >> rCentre_;
        periodicalOptionCoeffsDict.lookup("gradPAtflowDirection") >> gradPAtflowDirection_; 
        periodicalOptionCoeffsDict.lookup("Dh") >> Dh_; 

        Info << "\n   *The periodical option is selected as " << periodicalOption_ << endl;
        Info << "          The periodic angle is selected as " << thetaP_ << endl;             
        Info << "          The origin of the periodical axis is located at " << rCentre_ << endl;
        Info << "          The pressur gradient value induced in flow direction is " << gradPAtflowDirection_ << endl;
        Info << "          The length scale value used in ReDh is " << Dh_ << endl;

        //converting the thetaP to radian:
        thetaP_ = thetaP_ * mathematical::pi/180.0;
    }   
    else if( periodicalOption_ == "translationalPeriodicalInCurvedDuct")
    {
        periodicalOptionCoeffsDict.lookup("tPLength") >> tPLength_; 
        periodicalOptionCoeffsDict.lookup("rCentre") >> rCentre_;
        periodicalOptionCoeffsDict.lookup("gradPAtflowDirection") >> gradPAtflowDirection_; 
        periodicalOptionCoeffsDict.lookup("Dh") >> Dh_; 
 
        Info << "\n   *The periodical option is selected as " << periodicalOption_ << endl;
        Info << "          The periodic length is selected as " << tPLength_ << endl;
        Info << "          The origin of the periodical axis is located at " << rCentre_ << endl;
        Info << "          The pressur gradient value induced in flow direction is " << gradPAtflowDirection_ << endl;
        Info << "          The length scale value used in ReDh is " << Dh_ << endl;
    }
    else if( periodicalOption_ == "combinedRotationalPeriodicalInCurvedDuct" || periodicalOption_ == "combinedTranslationalPeriodicalInCurvedDuct" )
    {
        periodicalOptionCoeffsDict.lookup("thetaP") >> thetaP_; 
        periodicalOptionCoeffsDict.lookup("tPLength") >> tPLength_;
        periodicalOptionCoeffsDict.lookup("rCentre") >> rCentre_;
        periodicalOptionCoeffsDict.lookup("gradPAtflowDirection") >> gradPAtflowDirection_; 
        periodicalOptionCoeffsDict.lookup("Dh") >> Dh_; 

        Info << "\n   *The periodical option is selected as " << periodicalOption_ << endl;
        Info << "          The periodic angle is selected as " << thetaP_ << endl;             
        Info << "          The periodic length is selected as " << tPLength_ << endl;  
        Info << "          The origin of the periodical axis is located at " << rCentre_ << endl;
        Info << "          The pressure gradient value induced in flow direction is " << gradPAtflowDirection_ << endl;
        Info << "          The length scale value used in ReDh is " << Dh_ << endl;

        //converting the thetaP to radian:
        thetaP_ = thetaP_ * mathematical::pi/180.0;             
    }
    else if( periodicalOption_ == "straightDuctPeriodical")
    {
        periodicalOptionCoeffsDict.lookup("sDPLength") >> sDPLength_;
        periodicalOptionCoeffsDict.lookup("rCentre") >> rCentre_;
        periodicalOptionCoeffsDict.lookup("gradPAtflowDirection") >> gradPAtflowDirection_; 
        periodicalOptionCoeffsDict.lookup("Dh") >> Dh_; 

        Info << "\n   *The periodical option is selected as " << periodicalOption_ << endl;           
        Info << "          The periodic length  of straigh duct is selected as " << sDPLength_ << endl;  
        Info << "          The flow direction is parallel to the z axis" << endl;
        Info << "          The pressur gradient value induced in flow direction is " << gradPAtflowDirection_ << endl;
        Info << "          The length scale value used in ReDh is " << Dh_ << endl;             
    }    
    else if( periodicalOption_ == "none")
    {
        periodicalOptionCoeffsDict.lookup("gradPAtflowDirection") >> gradPAtflowDirection_; 
        periodicalOptionCoeffsDict.lookup("Dh") >> Dh_;     
        Info << "\n   *The periodical option is selected as " << periodicalOption_ << endl;
        Info << "          This means, there is no periodic B.C. in this problem." << endl;
        Info << "          but, we use a length scale, Dh for the ReDh calculation as: " << Dh_ << endl;
        Info << "          and the pressur gradient value induced in flow direction is: " << gradPAtflowDirection_ << endl;

    }

    // calc the important input parameters of the front mesh.
    lengthScaleOfTheMesh_ = calcLengthScaleOfTheMesh();
    if( frontMeshInputOption_ == "roughlyCalculated" || frontMeshInputOption_ == "finerCalculated" )
    {
        finerCalculatedMeshInputCoeffsDict.lookup("factorForMinEdge") >> factorForMinEdge_;
        finerCalculatedMeshInputCoeffsDict.lookup("factorForMaxEdge") >> factorForMaxEdge_;
        finerCalculatedMeshInputCoeffsDict.lookup("factorForMaxAspectRatio") >> factorForMaxAspectRatio_;

        minEdge_ = factorForMinEdge_ * lengthScaleOfTheMesh_;
        maxEdge_ = factorForMaxEdge_ * lengthScaleOfTheMesh_;
        maxAspectRatio_ = factorForMaxAspectRatio_ * 1.5;              

        Info << "\n   *All mesh parameters are calculated as: " << frontMeshInputOption_ << endl;
        Info << "       *The lengh scale of the mesh is " << lengthScaleOfTheMesh_ << endl;
        Info << "       *The minEdge parameter and it's factor for the front mesh is " << minEdge_  << ' ' << factorForMinEdge_ << endl;
        Info << "       *The maxEdge parameter and it's factor for the front mesh is " << maxEdge_  << ' ' << factorForMaxEdge_ << endl;
        Info << "       *The maxAspectRatio parameter and it's factor for the front mesh is " << maxAspectRatio_  << ' ' << factorForMaxAspectRatio_ << endl;

        if( frontMeshInputOption_ == "finerCalculated" )
        {
            finerCalculatedMeshInputCoeffsDict.lookup("howFinner") >> howFinner_;
            Info << "       *The howFinner parameter for finerCalculated option is " << howFinner_ << endl;
        }
    }     
    else if( frontMeshInputOption_ == "manually")
    {
        manuallyMeshInputCoeffsDict.lookup("howFinner") >> howFinner_; 
        manuallyMeshInputCoeffsDict.lookup("minEdge") >> minEdge_;
        manuallyMeshInputCoeffsDict.lookup("maxEdge") >> maxEdge_; 
        manuallyMeshInputCoeffsDict.lookup("maxAspectRatio") >> maxAspectRatio_; 

        Info << "\n   *All mesh parameters are selected as: " << frontMeshInputOption_ << endl;
        Info << "       *The calculated lengh scale of the mesh is " << lengthScaleOfTheMesh_ << endl;
        Info << "       *The value of howFinner showing how fine the front mesh should be is " << howFinner_ << endl;
        Info << "       *The minEdge parameter for the front mesh is " << minEdge_ << endl;
        Info << "       *The maxEdge parameter for the front mesh is " << maxEdge_ << endl;
        Info << "       *The maxAspectRatio parameter for the front mesh is " << maxAspectRatio_ << endl;             
    }


    // the important input parameter for smoothing and spreading the sources from front to Euelerian mesh.
    if( smoothingAndSpreadingOption_ == "calculated" )
    {   
        smoothingAndSpreadingCalculatedCoeffs.lookup("hFactor") >> hFactor_;
        h_ = hFactor_ * lengthScaleOfTheMesh_;     
        Info << "\n   *The parameter h for  the spreading and smoothing is selected as: " << smoothingAndSpreadingOption_ << endl;
        Info << "       *The lengh scale of the mesh is " << lengthScaleOfTheMesh_ << endl;
        Info << "       *The factor multyplying the h, hFactor is " << hFactor_ << endl;
        Info << "       *The parameter h is as " << h_ << endl;
    }     
    else if( smoothingAndSpreadingOption_ == "manually")
    {
        smoothingAndSpreadingManuallyCoeffs.lookup("h") >> h_;
        Info << "\n   *The parameter h for  the spreading and smoothing is selected as: " << smoothingAndSpreadingOption_ << endl;
        Info << "       *The lengh scale of the mesh is " << lengthScaleOfTheMesh_ << endl;
        Info << "       *The parameter h is selected as " << h_ << endl;          
    }

    // the important input parameter for constructing the indicators.
    if( constructingTheIndicatorsOption_ == "basedOnCPT" || constructingTheIndicatorsOption_ == "basedOnPoissonEquation" )
    {   
        constructingTheIndicatorsOptionCoeffs.lookup("numberOfIndicatorFiltering") >> numberOfIndicatorFiltering_;
        constructingTheIndicatorsOptionCoeffs.lookup("numberOfSearchToObtainTheExactCPForEachElement") >> numberOfSearchToObtainTheExactCPForEachElement_;
        constructingTheIndicatorsOptionCoeffs.lookup("filteringWeight") >> filteringWeight_;

        Info << "\n   *The method for constructing the indicators is selected as: " << constructingTheIndicatorsOption_ << endl;
        Info << "       *The number of indicator filtering is " << numberOfIndicatorFiltering_ << endl;
        Info << "       *The number of searching to obtain the exact closest point for each element is "
             << numberOfSearchToObtainTheExactCPForEachElement_ << endl;
        Info << "       *The filtering weight is selected as: " << filteringWeight_ << endl;
    }     
    else if ( constructingTheIndicatorsOption_ != "basedOnCPT" && constructingTheIndicatorsOption_ != "basedOnPoissonEquation" )
    {
        FatalErrorIn ("Foam::FrontDynamic<CloudType>::FrontDynamic ......")
            << "you should choose basedOnCPT or basedOnPoissonEquation for the constructingTheIndicatorsOption" << nl
            << abort(FatalError);
    }

    // the important input parameter for surface tension force distribution.
    if( surfaceTensionForceDistributionOption_ == "basedOnElementNeighbours" )
    {   
        Info << "\n   *The method for surface tension force distribution is selected as: " << surfaceTensionForceDistributionOption_ << endl;
    }
    else if( surfaceTensionForceDistributionOption_ == "basedOnPointNeighbours" )
    {   
        Info << "\n   *The method for surface tension force distribution is selected as: " << surfaceTensionForceDistributionOption_ << endl;
    }
    else
    {
        FatalErrorIn ("Foam::FrontDynamic<CloudType>::FrontDynamic ......")
            << "you should choose basedOnElementNeighbours or basedOnPointNeighbours for the surfaceTensionForceDistributionOption" << nl
            << abort(FatalError);
    }

    // the important input parameter for surface tension force distribution.
    if( undulationRemovalingOption_ == true )
    {   
        Info << "\n   *The method for undulation removaling operation is selected as: " << undulationRemovalingType_ << endl;
        Info << "       *The number of undulation removaling interval is " << undulationRemovalingInterval_ << endl;
    }   

    // the important input parameter for surface tension force distribution.
    if( pressureJumpAtTheInterfaceOption_ == true )
    {   
        Info << "\n   *The operation for reducting the spurious current will be done" << endl;
    } 

    flowDirectionField();
    //cacheFields();
    adjustingTheImportantVarsForIndicators();
    readInitialBubbles();
}

template<class CloudType>
Foam::FrontDynamic<CloudType>::FrontDynamic
(
    FrontDynamic<CloudType>& ft
)
:
    FrontAnalyzingModel<CloudType>(ft),
    ranGen_(ft.ranGen_),
    densityIndicator_(ft.densityIndicator_),
    viscosityIndicator_(ft.viscosityIndicator_),
    Indicator_(ft.Indicator_),
    cellsToRefine_(ft.cellsToRefine_)
{
}
// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //
template<class CloudType>
Foam::FrontDynamic<CloudType>::~FrontDynamic()
{}
// ************************************************************************* //
