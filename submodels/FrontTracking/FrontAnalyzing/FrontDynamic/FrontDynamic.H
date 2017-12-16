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

Class
    Foam::FrontDynamic

Description
    Collision model by 


\*---------------------------------------------------------------------------*/

#ifndef FrontDynamic_H
#define FrontDynamic_H

#include "FrontAnalyzingModel.H"
//Jafari added
#include "fileName.H"
#include "pointField.H"
#include "point.H"
#include "DLList.H"
#include "sphereToCell.H"
#include "cellSet.H"
#include "simpleMatrix.H"
#include "mathematicalConstants.H"//e,pi,twoPi,piByTwo
#include "fvCFD.H"

/*
//used for using the triSurface
#include "triSurface.H"
#include "PrimitivePatch.H"
#include "orientedSurface.H"
#include "triSurfaceSearch.H"
*/

#include "Random.H"
#include "IFstream.H"
#include "OFstream.H"
#include "IStringStream.H"//Jafari added
#include "OStringStream.H"//Jafari added
//#include "IOmanip.H"//Jafari added

#include <fstream> //Jafari added
#include <iostream> //Jafari added
#include <string>//Jafari added
#include <cmath>//Jafari added
#include <iomanip>//Jafari added
using namespace Foam::constant;//Jafari added
using std::fstream;//Jafari added
using std::ifstream;//Jafari added
using std::ofstream;//Jafari added
using std::iostream;//Jafari added
using std::ios;//Jafari added
using std::string;//Jafari added
using std::setw;//Jafari added


#include "UDictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                      Class FrontDynamic Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class FrontDynamic
:
    public FrontAnalyzingModel<CloudType>
{
protected:

    // Protected Data

        //- Convenience typedef to the cloud's parcel type
        typedef typename CloudType::parcelType parcelType;

        /*
        //Jafari added- to use triSurface
        typedef typename triSurface::elementInfo elTris;
        typedef typename triSurface::pointData ptTris;*/


        // Jafari added
        class pointData
        {
           public:
                  // The index of the current point. 
                  label currentIndex;

                  // The list of connected points to the current point as a DynamicList<pointData>. 
                  //DynamicList<pointData> connectedPoints;

                  DynamicList<point> connectedPoints;
                  DynamicList<label> connectedPointsIndex;

                  // positions of the current point as a vector. 
                  vector currentPoint;
                  vector posInDomain;

                  // this is used when an element must be eliminated from DynamicList of elementData.
                  bool keepPoint;

                  label elOwner;

                  pointData()
                  {
                      currentIndex = -1;
                      currentPoint = vector::zero;
                      posInDomain = vector::zero;
                      keepPoint = true;
                      elOwner = -1;
                  };
        };

        // Jafari added
        class elementInfo
        {
           public:
                  // The index of the current element. 
                  label currentIndex;

                  // indexes of the elements which are neighbours of the current one in counter clockwise order. 
                  labelList elementIndex;

                  // indexes of the points which are ends of the current element in counter clockwise order. 
                  labelList pointIndex;

                  // position of the points which are ends of the current element in counter clockwise order. 
                  List<point> points;

                  // this is used when an element must be eliminated from DynamicList of elementData.
                  bool keepElement;

                  // required data 
                  vector centrePosInDomain;
                  vector elementSurfaceArea;
                  vector velocity;
                  scalar pressure;
                  symmTensor shearStress;
                  vector averageSurfaceTension; 

                  elementInfo()
                  {
                      elementIndex.setSize(3);
                      pointIndex.setSize(3);
                      points.setSize(3);
                      currentIndex = -1;
                      for( int i=0; i<=2; i++)
                      {
                          elementIndex[i] = -1;
                          pointIndex[i] = -1;
                          points[i] = vector::zero;
                      }
                      keepElement = true;
                      centrePosInDomain = vector::zero;
                      elementSurfaceArea = vector::zero;
                      velocity = vector::zero;
                      pressure = 0.0;
                      shearStress = symmTensor::zero;
                      averageSurfaceTension = vector::zero;
                  };

        };

        // Jafari added
        class CPTData // closest point transform data
        {
           public:
                  // list of density and viscosity induced by each bubble. 
                  List<scalar> inducedDensity;
                  List<scalar> inducedViscosity;
                  CPTData()
                  {
                  };
        };

        // Jafari added
        class bubbleData
        {
           public: 
                  // The pointData List for the current bubble. 
                  DynamicList<pointData> pL;

                  // The elementInfo List for the current bubble. 
                  DynamicList<elementInfo> eL;

                  // The pressure force exerted on the current bubble. 
                  vector pressureForce;

                  // The shear stress (friction force) force exerted on the current bubble. 
                  vector shearStressForce;

                  // The pressure drag force value exerted on the current bubble. 
                  scalar pressureDragForce;

                  // The friction drag force value exerted on the current bubble. 
                  scalar frictionDragForce;

                  // The pressure lift force value exerted on the right side of the current bubble. 
                  scalar pressureLiftForceI;

                  // The pressure lift force value exerted on the last side of the current bubble. 
                  scalar pressureLiftForceII;

                  // The friction lift force exerted on the right side of the current bubble. 
                  scalar frictionLiftForceI;

                  // The friction lift force exerted on the last side of the current bubble. 
                  scalar frictionLiftForceII;

                  // The total drag force exerted on the current bubble. 
                  scalar totalDragForce;

                  // The total lift force value exerted on the right side of the current bubble. 
                  scalar totalLiftForceI;

                  // The total lift force value exerted on the last side of the current bubble.  
                  scalar totalLiftForceII;

                  // The buoyant force exerted on the current bubble. 
                  vector buoyantForce;

                  // The gravity force exerted on the current bubble. 
                  vector gravityForce;

                  // The accelaration of the current bubble. 
                  vector accelaration;

                  // The velocity of the current bubble. 
                  vector velocity;

                  // The velocity of the current bubble at old time. 
                  vector velocityOld;

                  // The centre of the current bubble. 
                  vector centre;

                  // The centre of the current bubble at old time. 
                  vector centreOld;

                  // The volume of the current bubble. 
                  scalar volume;

                  // The surface area of the current bubble. 
                  scalar surfaceArea;

                  // The current bubble mass. 
                  scalar mass;

                  // The bubble density. 
                  scalar density;

                  // The bubble density ratio. 
                  scalar outerFluidDensity;

                  // The bubble viscosity. 
                  scalar viscosity;

                  // The bubble viscosity ratio. 
                  scalar outerFluidViscosity;

                  // The bubble drag coefficient. 
                  scalar dragCoeff;

                  // The bubble lift coefficient of right side of motion direction. 
                  scalar LiftCoeffI;

                  // The bubble lift coefficient of last side direction. 
                  scalar LiftCoeffII;

                  // The index of the current bubble. 
                  label currentIndex;

                  // The surface tension Coeff of the current bubble. 
                  scalar surfaceTensionCoeff;

                  // The Morton number of the current bubble. 
                  scalar MortonNumber;

                  // The Eotvos number of the current bubble. 
                  scalar EotvosNumber;

                  // The equalized diameter of the current bubble when it has the same volume with a sphere. 
                  scalar sphereBubbleDiameter;

                  // The maximum diameter, each diameter component show the maximum value of bubble in that direction. 
                  List<scalar> diameter;

                  // The sphericity showing how the bubble shape diffes from a sphere with the same volume. 
                  scalar sphericity;

                  vector pressureJumpAtTheInterface;

                  // this is used when an bubble must 
                  // be eliminated from DynamicList of bubbleData.
                  bool keepBubble; 

                  bubbleData()
                  {
                      // The pointData List for the current bubble. 
                      pL.clearStorage();
                      eL.clearStorage();
                      pressureForce = vector::zero;
                      shearStressForce = vector::zero;
                      pressureDragForce = 0.0;
                      frictionDragForce = 0.0;
                      pressureLiftForceI = 0.0;
                      pressureLiftForceII = 0.0;
                      frictionLiftForceI = 0.0;
                      frictionLiftForceII = 0.0;
                      totalDragForce = 0.0;
                      totalLiftForceI = 0.0;
                      totalLiftForceII = 0.0;
                      buoyantForce = vector::zero;
                      gravityForce = vector::zero;
                      accelaration = vector::zero;
                      velocity = vector::zero;
                      velocityOld = vector::zero;
                      centre = vector::zero;
                      centreOld = vector::zero;
                      volume = 0.0;
                      surfaceArea = 0.0;
                      mass = 0.0;
                      density = 0.0;
                      outerFluidDensity = 0.0;
                      viscosity = 0.0;
                      outerFluidViscosity = 0.0;
                      dragCoeff = 0.0;
                      LiftCoeffI = 0.0;
                      LiftCoeffII = 0.0;
                      currentIndex = 0.0;
                      surfaceTensionCoeff = 0.0;
                      MortonNumber = 0.0;
                      EotvosNumber = 0.0;
                      sphereBubbleDiameter = 0.0;
                      diameter.setSize(3);
                      sphericity = 0.0;
                      pressureJumpAtTheInterface = vector::zero;
                      keepBubble = true;
                  };
        };
       
    // The bubbleData list of the bubbleData class.

        DynamicList<bubbleData> bDataL_; 

	word initialMakingFileOfTheBubbles_;

        //- Random number generator
        Random ranGen_;

        scalar timeTobeSteady_ = 0;

        scalar timeOld_ = 0.0; // container of old time

    // The periodical options for curved duct. note: At theta =0, the r( or polar) axis
    // is the very x-axis and y axis is alinged to flow direction

        //- the centre of periodical axis ( z-axis can be rotational and/or translational periodical axis)
        vector rCentre_ = vector::zero;

        //- the periodical angle about rotational axis
        scalar thetaP_ = 0.0;

        //- the translational periodic length in the direction of z-axis
        scalar tPLength_ = 0.0;

        //- the periodic length of straight duct in the direction of z-axis
        scalar sDPLength_ = 0.0;

        //- periodical option
        word periodicalOption_;// rotationalPeriodicalInCurvedDuct/translationalPeriodicalInCurvedDuct/
                               // combinedRotationalPeriodicalInCurvedDuct/combinedTranslationalPeriodicalInCurvedDuct/none/straightDuctPeriodical

        vector flowD_ = vector::zero; // flow direction.
        vector normalToflowDI_ = vector::zero; // the first (radial direction) normal to flow direction.
        vector normalToflowDII_ = vector::zero; // the second (z direction) normal to flow direction.
        DynamicList<vector> flowDList_;// list of all flow direction at each cell.


        //scalar EotvosNumber_ = 0.0;
        //scalar MortonNumber_ = 0.0;
        scalar ReDh_ = 0.0;
        scalar Dh_ = 0.0; //length used in ReDh
        scalar surfaceTensionCoeff_ = 0.0;
        scalar magFlowVelocity_ = 0.0;
        scalar gradPAtflowDirection_ = 0.0;

    // The important input parameters of the front mesh.

        //- periodical option
        word frontMeshInputOption_;// roughlyCalculated/finerCalculated/manually

        scalar lengthScaleOfTheMesh_;

        scalar minEdge_;

        scalar maxEdge_;

        scalar maxAspectRatio_;

        labelList locInPointDataList_;

        labelList locInElementInfoList_;

        scalar howFinner_; // in the range 0.05 (high finner)- 0.3 (low finner)

        scalar factorForMinEdge_; // in the range 0.2 (high finner)- 0.3 (low finner)
        scalar factorForMaxEdge_; // in the range 1.0 (high finner)- 2.0 (low finner)
        scalar factorForMaxAspectRatio_; // in the range 1 (high finner)- 1.33 (low finner)

    // undulation removaling operation switch

       bool  undulationRemovalingOption_ = false;
       label undulationRemovalingInterval_ = 0;
       word  undulationRemovalingType_; // TSUR3D/volumeConservingSmoothingIII/volumeConservingSmoothingIV
       label URIntervalCounter_ = 0;
       scalar initialBubbleVolume_;

    // constructing the indicator functions

       List<label> p_; 
       List<label> q_; 
       List<label> r_;
       //simpleMatrix<scalar> mComponent_;
       List<scalar> beta_;
       scalar h_;
       scalar hFactor_;     
       word smoothingAndSpreadingOption_; // manually/calculated
       volScalarField densityIndicator_;
       volScalarField viscosityIndicator_;
       volScalarField Indicator_;
       volScalarField cellsToRefine_;

       List<CPTData> CPTList_;
       List<DynamicList<label> > allCellsInEeachMasket_;
       List<DynamicList<label> > allCellsNearEeachFront_;
       List<scalar> lengthScaleNearTheFront_;
       word constructingTheIndicatorsOption_; // manually/calculated
       label numberOfIndicatorFiltering_ = 0;
       label numberOfSearchToObtainTheExactCPForEachElement_ =0;
       scalar filteringWeight_ = 0.0;


    // for surface tension force distribution
       word surfaceTensionForceDistributionOption_; // basedOnElementNeighbours/basedOnPointNeighbours

    // for pressure jump at the interface to reduce spurious currents
       bool  pressureJumpAtTheInterfaceOption_ = false;

    // operations for phisical margins to density and viscosity
       bool  twoFluidFlow_ = false;

    // Protected Data

	//- p_rgh
        const volScalarField* p_rghPtr_;

	//- p
        const volScalarField* pPtr_;

	//- U
        const volVectorField* UPtr_;

    // search in box activation for reduction the computational cost
       bool  searchInBox_ = true;


    // Protected Member Functions

        //- Main collision routine
        virtual void dynamic();

        //- Jafari added
        void readInitialBubbles();
    
public:

    //- Runtime type information
    TypeName("FrontDynamic");


    // Constructors

        //- Construct from dictionary
        FrontDynamic
        (
            const dictionary& dict,
            CloudType& cloud,
            const word& modelName = typeName
        );

        //- Construct copy
        FrontDynamic(FrontDynamic<CloudType>& ft);

        //- Construct and return a clone
        virtual autoPtr<FrontAnalyzingModel<CloudType> > clone() //const
        {
            return autoPtr<FrontAnalyzingModel<CloudType> >
            (
                new FrontDynamic<CloudType>(*this)
            );
        }

    //- Destructor
        virtual ~FrontDynamic();

    // Member Functions

        inline void cloudCleaning();

        inline void mapFrontToParcel();

	    inline void calculatingThePosInDomain();

	    inline Foam::vector calcThePositionInDomain( point pos);

	    inline void calculatingThePosInDomainOfPoints();

	    inline void calculatingTheCentrePosInDomainOfElements();

        inline void mapParcelToFront();

	template <class T>
        inline Foam::label quickSelect(const DynamicList<T>& S, const label k, label startI, label endI);

        inline void makingBubblesAtTheZeroTime();

        inline void findNighboures(DynamicList<pointData>& pL, DynamicList<elementInfo>& eL);

        inline void pointToElementMapping(DynamicList<pointData>& pL, DynamicList<elementInfo>& eL);

        inline void ptConnectedPoints( DynamicList<pointData>& pL, DynamicList<elementInfo>& eL);

        inline void findPointsOrder( label f, label& cNI, label& cNII, label& cNIII );

        inline void printInitialFronts(const fileName& outfileName, DynamicList<pointData> pL, DynamicList<elementInfo> eL);

        inline void startPtIndexFromZero(DynamicList<pointData>& pL, DynamicList<elementInfo>& eL);

        inline bool isElementSmall(elementInfo el, scalar minEdge, scalar maxEdge,
                                   scalar maxAspectRatio, scalar& s1, scalar& s2,
                                   scalar& s3, label& n1, label& n2, label& n3);

	    inline bool isElementRefined(elementInfo el, scalar minEdge, scalar maxEdge, scalar maxAspectRatio,
                                     scalar& s1, scalar& s2, scalar& s3, label& n1, label& n2, label& n3);

	    inline void frontCoarsening(DynamicList<pointData>& ptDataList, DynamicList<elementInfo>& elInfoList,
                                                                   scalar minEdge, scalar maxEdge, scalar maxAspectRatio);

	    inline void frontRefining( DynamicList<pointData>& ptDataList, DynamicList<elementInfo>& elInfoList,
                                   scalar minEdge, scalar maxEdge, scalar maxAspectRatio);

	    inline void TSUR3D( label repeatNum, DynamicList<pointData>& ptDataList, DynamicList<elementInfo>& elInfoList);

	    inline void volumeConservingSmoothingIII( label repeatNum, DynamicList<pointData>& ptDataList, DynamicList<elementInfo>& elInfoList);

	    inline void volumeConservingSmoothingIV( label repeatNum, DynamicList<pointData>& ptDataList, DynamicList<elementInfo>& elInfoList);

	    inline Foam::vector smoothTo( pointData currentPt);

    	inline Foam::scalar calcBubbleVolume( DynamicList<elementInfo> eL);

	   inline Foam::scalar calcBubbleSurface( DynamicList<elementInfo> eL);

	   inline Foam::scalar PeskinD(scalar r);

	   inline Foam::scalar windowFunc(vector del, scalar h);

	   inline Foam::scalar windowFuncBar(vector del, scalar h, List<scalar> beta);

	   inline Foam::scalar mPQR(vector del, scalar h, label p, label q, label r);

	   inline Foam::vector calcTheSurfaceAreaAtpoint( pointData ptD);

       inline Foam::vector calcSurfaceTensionForceAtpoint( pointData ptD);

	   inline Foam::vector calcElementSurfaceArea( elementInfo elD);

	static string getLineNoComment(IFstream& is)
	{
	    string line;
	    do
	    {
		is.getLine(line);
	    }
	    while ((line.empty() || line[0] == '#') && is.good());
	    return line;
	}

        inline void printFrontsForAllTimes();

        inline void writeDataOfTheFrontsAtTheOutputTime();

        inline void writeDataOfTheFrontsAtTheOutputTimeIntimeName();//JJ

        inline void readDataOfTheFrontsAtTheCurrentTime();

        inline void readDataOfTheFrontsAtTheCurrentTimeIntimeName();//JJ

        inline Foam::scalar calcLengthScaleOfTheMesh();
   
        inline void processingThefrontMeshes();

        inline void undulationRemovalingOption();

        inline void adjustingTheMeshDataAfterChanges(DynamicList<pointData>& pL, DynamicList<elementInfo>& eL);

        inline void adjustingTheImportantVarsForIndicators();

        inline void constructingTheIndicatorsBasedOnCPT();

        inline void calculatingTheExactClosestPointForTheCurrentElement(elementInfo elInfo, vector& minDistance);

        inline void communicationsOfTheFrontAndEulerianGrid();

        inline void searchAllCellsInEeachMasket();

        inline void searchAllCellsNearEeachFront();

        inline void surfaceTensionForceDistribution();

        inline void cacheFields();

        inline Foam::tmp<Foam::volScalarField> p_rgh() const;

        inline Foam::tmp<Foam::volScalarField> p() const;

        inline Foam::tmp<Foam::volScalarField> rho() const;

        inline Foam::tmp<Foam::volScalarField> mu() const;

        inline Foam::tmp<Foam::volVectorField> U() const;

        inline Foam::tmp<Foam::volScalarField> densityIndicator();

        inline Foam::tmp<Foam::volScalarField> viscosityIndicator();

        inline void bubblePostProcessing();

        inline void writeBubblePositions();

        inline void writeBubbleVelocityAndAccelaration();

        inline void writeBubbleSphericity();

        inline void writeBubbleDragCoefs();

        inline void writeBubbleLiftCoefs();

        inline void writeBubbleDiameterComponents();

        inline void calcBubbleNonDimensionalNumbers();

        inline void outputTimePostProcessing();

        inline void flowDirectionField();

        inline void printingAllImportantRunTimeData();

        inline void initialAllocatingTheCPTList();

        inline Foam::scalar LiuHeviside( scalar outerFluidProperty, scalar property, vector minDistance, vector elementSurfaceArea, scalar gama);

        inline void setLengthScalesNearTheFront();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

 #include "functionsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "FrontDynamic.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
