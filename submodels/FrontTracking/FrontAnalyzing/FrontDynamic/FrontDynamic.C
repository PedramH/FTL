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

//--------------------------------------------------------------------------
//-----------------------------cloudCleaning--------------------------------
//--------------------------------------------------------------------------
template<class CloudType>
inline void Foam::FrontDynamic<CloudType>::cloudCleaning()
{
   forAllIter(typename CloudType, this->owner(), iter1)
   {
       parcelType& p1 = iter1();
       this->owner().deleteParticle( p1);
   }
   /*
   forAllIter(typename CloudType, this->owner(), iter1)
   {
       //parcelType& p1 = iter1();
       this->owner().deleteParticle( iter1());
   }
   */
    Info << "\n---------> At the time: " << this->owner().db().time().value()
         <<" , Now all front points is deleting ....... " << endl;
}


//--------------------------------------------------------------------------
//-----------------------------mapFrontToParcel--------------------------------
//--------------------------------------------------------------------------
template<class CloudType>
inline void Foam::FrontDynamic<CloudType>::mapFrontToParcel()
{
    Info << "\n---------> At the time: " << this->owner().db().time().value()
         <<" , Now all front points is mapping to the parcels ....... " << endl;

    const fvMesh& mesh = this->owner().mesh();

    // insert correctly the particles................
    if (Pstream::parRun())
    {//2

        forAll(bDataL_,bDI)
        {//10
            DynamicList<pointData> tempPL = bDataL_[bDI].pL;

	    // Determine the front point position and owner cell,
	    // tetFace and tetPt
	    label cellI = -1;
	    label tetFaceI = -1;
	    label tetPtI = -1;
	    vector pos = vector::zero;

	    label posInList = -1;

	    forAll(tempPL, pointI)
	    {//3
		pointData ptD = tempPL[pointI];
		pos = ptD.posInDomain;
		this->owner().mesh().findCellFacePt
		(
		    pos,
		    cellI,
		    tetFaceI,
		    tetPtI
		);

	        // gather cellI from all proc
		labelList LcellI(Pstream::nProcs());
		LcellI[Pstream::myProcNo()] = cellI;
		Pstream::gatherList(LcellI);
    		Pstream::scatterList(LcellI);

	        // gather tetFaceI from all proc
    		labelList LtetFaceI(Pstream::nProcs());
    		LtetFaceI[Pstream::myProcNo()] = tetFaceI;
    		Pstream::gatherList(LtetFaceI);
    		Pstream::scatterList(LtetFaceI);

	        // gather tetPtI from all proc
    		labelList LtetPtI(Pstream::nProcs());
    		LtetPtI[Pstream::myProcNo()] = tetPtI;
    		Pstream::gatherList(LtetPtI);
    		Pstream::scatterList(LtetPtI);

	    	forAll(LcellI,LcI)
	   	{
	     	    if( LcellI[LcI] > -1 && LtetFaceI[LcI] > -1 )
	            {
	                posInList = LcI;
	                break;
	            }
	        }

		if (  Pstream::myProcNo() == posInList )
		{
	             parcelType* pPtr = new parcelType(mesh, pos, LcellI[posInList], LtetFaceI[posInList], LtetPtI[posInList]);
	             //Check/set new parcel thermo properties
	             this->owner().setParcelThermoProperties(*pPtr, 0.0);
	             //Apply correction to velocity for 2-D cases
	             meshTools::constrainDirection(this->owner().mesh(), this->owner().mesh().solutionD(), pPtr->U());
	             //Check/set new parcel injection properties
	             this->owner().checkParcelProperties(*pPtr, 0.0, false);
	             //Apply correction to velocity for 2-D cases
	             meshTools::constrainDirection
		     (
			 this->owner().mesh(),
			 this->owner().mesh().solutionD(),
			 pPtr->U()
		     );
		     pPtr->currentIndex() = ptD.currentIndex;
   	             pPtr->bubbleIndex() = bDI;
		     pPtr->shadowPos() = ptD.currentPoint;
		     this->owner().addParticle(pPtr);
		 }
	    }//3

        }//10
    }//2

    if (!Pstream::parRun())
    {//5
        forAll(bDataL_,bDI)
        {//10
	    DynamicList<pointData> tempPL = bDataL_[bDI].pL;

	    // Determine the front point position and owner cell,
	    // tetFace and tetPt
	    label cellI = -1;
	    label tetFaceI = -1;
	    label tetPtI = -1;
	    vector pos = vector::zero;

	    forAll(tempPL, pointI)
	    {//4
		pointData ptD = tempPL[pointI];
		pos = ptD.posInDomain;
		this->owner().mesh().findCellFacePt
		(
		    pos,
		    cellI,
		    tetFaceI,
		    tetPtI
		);

		if ( cellI > -1)
		{
                     //Info << " AF $$$$$$$$$$ cellI/tetFaceI/tetPtI " << cellI << ' ' << tetFaceI << ' ' << tetPtI << endl;
	            parcelType* pPtr = new parcelType(mesh, pos, cellI, tetFaceI, tetPtI);
	    	    //Check/set new parcel thermo properties
	            this->owner().setParcelThermoProperties(*pPtr, 0.0);
	            //Apply correction to velocity for 2-D cases
	            meshTools::constrainDirection(this->owner().mesh(), this->owner().mesh().solutionD(), pPtr->U());
		    // Check/set new parcel injection properties
		    this->owner().checkParcelProperties(*pPtr, 0.0, false);
		    // Apply correction to velocity for 2-D cases
		    meshTools::constrainDirection
		    (
		        this->owner().mesh(),
			this->owner().mesh().solutionD(),
			pPtr->U()
		    );
		    pPtr->currentIndex() = ptD.currentIndex;
   	            pPtr->bubbleIndex() = bDI;
		    pPtr->shadowPos() = ptD.currentPoint;
		    this->owner().addParticle(pPtr);
		 }
	    }//4
        }//10
    }//5
    // insert correctly the particles................
}
//--------------------------------------------------------------------------
//--------------------------------calculatingThePosInDomain-----------------
//--------------------------------------------------------------------------
template<class CloudType>
inline void Foam::FrontDynamic<CloudType>::calculatingThePosInDomain()
{
    calculatingThePosInDomainOfPoints();
    calculatingTheCentrePosInDomainOfElements();
}

//note: func to calculate the position in domain.
template<class CloudType>
inline Foam::vector Foam::FrontDynamic<CloudType>::calcThePositionInDomain( point pos)
{
    vector posInDomain;

    if( periodicalOption_ == "rotationalPeriodicalInCurvedDuct")
    {//0
	rCentre_.z() = pos.z();
        vector rp = pos - rCentre_;
	scalar ptTheta = acos( rp.x()/mag(rp));
	scalar ptThetaInDomain = 0.0;
	if( rp.y() >= 0 )
	{
	    ptThetaInDomain = ptTheta - thetaP_ * floor(ptTheta/thetaP_);
	}
	else if( rp.y() < 0 )
	{
	    ptTheta = mathematical::twoPi - ptTheta;
	    ptThetaInDomain = ptTheta - thetaP_ * floor(ptTheta/thetaP_);
	}
	// calculating the posInDoamain
	posInDomain.x() = mag(rp) * cos(ptThetaInDomain);
	posInDomain.y() = mag(rp) * sin(ptThetaInDomain);
	posInDomain.z() = pos.z();
        return posInDomain;
    }//0
    else if( periodicalOption_ == "translationalPeriodicalInCurvedDuct")
    {
	rCentre_.z() = pos.z();
        vector rp = pos - rCentre_;
	// calculating the posInDoamain
	posInDomain = pos;
	posInDomain.z() = rp.z() - tPLength_ * floor(rp.z()/tPLength_);
        return posInDomain;
    }
    else if( periodicalOption_ == "combinedRotationalPeriodicalInCurvedDuct" || periodicalOption_ == "combinedTranslationalPeriodicalInCurvedDuct" )
    {
	rCentre_.z() = pos.z();
        vector rp = pos - rCentre_;
	scalar ptTheta = acos( rp.x()/mag(rp));
	scalar ptThetaInDomain = 0.0;
	if( rp.y() >= 0 )
	{
	    ptThetaInDomain = ptTheta - thetaP_ * floor(ptTheta/thetaP_);
	}
	else if( rp.y() < 0 )
	{
	    ptTheta = mathematical::twoPi - ptTheta;
	    ptThetaInDomain = ptTheta - thetaP_ * floor(ptTheta/thetaP_);
	}
	// calculating the posInDoamain
	posInDomain.x() = mag(rp) * cos(ptThetaInDomain);
	posInDomain.y() = mag(rp) * sin(ptThetaInDomain);
        // this is done for the translationalPeriodical part
	posInDomain.z() = rp.z() - tPLength_ * floor(rp.z()/tPLength_);
        return posInDomain;
    }
    else if( periodicalOption_ == "straightDuctPeriodical")
    {
	// calculating the posInDoamain
	posInDomain = pos;
	posInDomain.z() = pos.z() - sDPLength_ * floor(pos.z()/sDPLength_);
        return posInDomain;
    }
    else if( periodicalOption_ == "none")
    {
	// calculating the posInDoamain
	posInDomain = pos;
        return posInDomain;
    }
    else
    {
        FatalErrorIn ("inline Foam::vector Foam::FrontDynamic<CloudType>::calcThePositionInDomain( point pos)")
            << "you should choose periodicalOption correctly" << nl
            << abort(FatalError);
    }
}

//note: func to calculate the posInDomain of all points.
template<class CloudType>
inline void Foam::FrontDynamic<CloudType>::calculatingThePosInDomainOfPoints()
{
    Info <<"\n---------> At the time: " << this->owner().db().time().value()
         <<" , Now all front point locations is mapping to the domain ....... " << endl;
    forAll(bDataL_,bDI)
    {//1
        DynamicList<pointData>& pL = bDataL_[bDI].pL;
	forAll(pL,pI)
	{//2
	    point pos = pL[pI].currentPoint;
            pL[pI].posInDomain = calcThePositionInDomain(pos);
        }//2
    }//1
}

//note: func to calculate the posInDomain of all points.
template<class CloudType>
inline void Foam::FrontDynamic<CloudType>::calculatingTheCentrePosInDomainOfElements()
{
    Info << "\n---------> At the time: " << this->owner().db().time().value()
         <<" , Now all front element centres is mapping to the domain ....... " << endl;

    forAll(bDataL_,bDI)
    {//1
        DynamicList<pointData>& pL = bDataL_[bDI].pL;
	DynamicList<elementInfo>& eL = bDataL_[bDI].eL;
	forAll(eL,eI)
	{//2
	    label p1I = eL[eI].pointIndex[0];
	    label p2I = eL[eI].pointIndex[1];
	    label p3I = eL[eI].pointIndex[2];
	    point pos = (pL[p1I].currentPoint + pL[p2I].currentPoint + pL[p3I].currentPoint)/3.0;
            eL[eI].centrePosInDomain = calcThePositionInDomain(pos);
        }//2
    }//1
}
//--------------------------------------------------------------------------
//--------------------------------mapParcelToFront--------------------------
//--------------------------------------------------------------------------
template<class CloudType>
inline void Foam::FrontDynamic<CloudType>::mapParcelToFront()
{
    Info << "\n---------> At the time: " << this->owner().db().time().value()
         <<" , Now all parcels is mapping to the front ....... " << endl;

    if (Pstream::parRun())
    {//2
	DynamicList<label>  myTranferLabelList;
	DynamicList<label>  myTranferBLabelList;
	DynamicList<vector> myTranferPointList;

	forAllIter(typename CloudType, this->owner(), iter1)
	{
	    parcelType& p1 = iter1();
	    myTranferLabelList.append( p1.currentIndex());
	    myTranferBLabelList.append( p1.bubbleIndex());
	    myTranferPointList.append( p1.shadowPos());
	}

	DynamicList<label>  totalTranferLabelList;
	DynamicList<label>  totalTranferBLabelList;
	DynamicList<vector> totalTranferVectorList;

	// gather/scatter for myTranferLabelList
	List<DynamicList<label> > allLL(Pstream::nProcs());
	allLL[Pstream::myProcNo()] = myTranferLabelList;
	Pstream::gatherList(allLL);
	Pstream::scatterList(allLL);
	// gather/scatter for myTranferBLabelList
	List<DynamicList<label> > allBLL(Pstream::nProcs());
	allBLL[Pstream::myProcNo()] = myTranferBLabelList;
	Pstream::gatherList(allBLL);
	Pstream::scatterList(allBLL);
	// gather/scatter for myTranferPointList
	List<DynamicList<vector> > allVL(Pstream::nProcs());
	allVL[Pstream::myProcNo()] = myTranferPointList;
	Pstream::gatherList(allVL);
	Pstream::scatterList(allVL);

	// preparing the totalTranferLabelList & totalTranferVectorList
	forAll(allBLL,BLI)
	{
	    DynamicList<label> tempLL = allLL[BLI];
	    DynamicList<label> tempBLL = allBLL[BLI];
	    DynamicList<vector> tempVL = allVL[BLI];

	    forAll(tempBLL,tBLI)
	    {
		totalTranferLabelList.append(tempLL[tBLI]);
		totalTranferBLabelList.append(tempBLL[tBLI]);
		totalTranferVectorList.append(tempVL[tBLI]);
	    }
	}

	// updating the pL
	forAll(totalTranferBLabelList,tTBLI)
	{
            label bubbleDI = totalTranferBLabelList[tTBLI];
            DynamicList<pointData>& pL = bDataL_[bubbleDI].pL;

	    //label posInList = quickSelect(pL, totalTranferLabelList[tTBLI], 0, pL.size()-1);
	    //pointData& ptD = pL[posInList];
	    pointData& ptD = pL[totalTranferLabelList[tTBLI]];
	    ptD.currentPoint = totalTranferVectorList[tTBLI];
	}

	// updating the eL
	forAll(bDataL_,bDI)
	{
            DynamicList<pointData>& pL = bDataL_[bDI].pL;
            DynamicList<elementInfo>& eL = bDataL_[bDI].eL;
	    pointToElementMapping(pL,eL);
	}

	totalTranferLabelList.clearStorage();
	totalTranferVectorList.clearStorage();
	myTranferLabelList.clearStorage();
	myTranferPointList.clearStorage();
    }//2

    if (!Pstream::parRun())
    {//3
	forAllIter(typename CloudType, this->owner(), iter1)
	{
	    parcelType& p1 = iter1();

            label bubbleDI = p1.bubbleIndex();
            label pointDI = p1.currentIndex();

            DynamicList<pointData>& pL = bDataL_[bubbleDI].pL;
	    //label posInList = quickSelect(pL, pointDI, 0, pL.size()-1);
	    //pointData& ptD = pL[posInList];
	    pointData& ptD = pL[pointDI];

	    ptD.currentPoint = p1.shadowPos();
	}

	// updating the eL
	forAll(bDataL_,bDI)
	{
            DynamicList<pointData>& pL = bDataL_[bDI].pL;
            DynamicList<elementInfo>& eL = bDataL_[bDI].eL;
	    pointToElementMapping(pL,eL);
	}
    }//3
}

//--------------------------------------------------------------------------
//-----------------------------quickSelect----------------------------------
//--------------------------------------------------------------------------
// quichSelect func is used to find positon of any member in
// DynamicList<elementData> with currentIndex.
template<class CloudType>
template <class T>
inline Foam::label Foam::FrontDynamic<CloudType>::quickSelect(const DynamicList<T>& S, const label k, label startI, label endI)
{
    scalar rand = ranGen_.scalar01();
    label pivot = startI + label( rand * ( endI - startI) );

    if( S[pivot].currentIndex < k )
    {
        startI = pivot + 1;
        return quickSelect(S,k, startI, endI);
    }
    else if ( S[pivot].currentIndex > k )
    {
        endI   = pivot -1;
        return quickSelect(S,k, startI, endI);
    }
    else
    {
        return pivot;
    }
}

//--------------------------------------------------------------------------
//-----------------------------makingBubblesAtTheZeroTime--------------------------------------
//--------------------------------------------------------------------------
//note: func to read the all bubble front mesh.
template<class CloudType>
inline void Foam::FrontDynamic<CloudType>::makingBubblesAtTheZeroTime()
{
    fileName infileName = initialMakingFileOfTheBubbles_;

    // This is the directory keeping the bubble data.
    fileName FTResult = "FTResult";

    fileName InfileNamePath =  this->owner().db().time().rootPath()/this->owner().db().time().globalCaseName()/FTResult/infileName;

    IFstream inFile(InfileNamePath);

    label bubbleNum;
    scalar radin;
    label nps;
    scalar xc, yc, zc, e;// bubble centre
    scalar outerFluidDensity, density, outerFluidViscosity, viscosity;// bubble properties
    scalar sigma;// surface tension coeff
    scalar Eotvos;// Eotvos number for the current bubble
    scalar Morton;// Morton number for the current bubble

    pointData ptI; //tempoaray pointData var
    elementInfo elI; //tempoaray elementInfo var

    string line = getLineNoComment(inFile);
    IStringStream lineStream(line);
    lineStream >> bubbleNum;

    Info << "\n   *Bubble construction process ....." << endl;
    Info << "\n       The total number of bubbles is " << bubbleNum << endl;

    List<fileName> printAllbubbles(bubbleNum);
    forAll(printAllbubbles,pAB)
    {
        // note, this way is for intrinsic C++ utility
	//std::ostringstream oss;
	//oss << "b" << pAB << ".plt";
	//printAllbubbles[pAB] = oss.str();

        // note, this way is special for OpenFOAM using the OStringStream object
        OStringStream oss;
	oss << "bubble" << pAB << "_AtZeroTime" << ".plt";
        fileName fName = oss.str();
	printAllbubbles[pAB] =  this->owner().db().time().rootPath()/this->owner().db().time().globalCaseName()/FTResult/fName;
    }

    for(int bI = 0; bI < bubbleNum ; bI++)
    {//0
        bubbleData bData;
        DynamicList<pointData>& pL = bData.pL;
        DynamicList<elementInfo>& eL = bData.eL;

        pL.reserve( 100000);//JJ
        eL.reserve( 200000);//JJ

        string line = getLineNoComment(inFile);
        {
           IStringStream lineStream(line);
           lineStream >> radin >> nps >> xc >> yc >> zc >> e
                      >> outerFluidDensity >> density
                      >> outerFluidViscosity >> viscosity
                      >> sigma >> Eotvos >> Morton;
        }

        bData.outerFluidDensity = outerFluidDensity;
        bData.density = density;
        bData.outerFluidViscosity = outerFluidViscosity;
        bData.viscosity = viscosity;
        bData.sphereBubbleDiameter = 2.0 * radin;
        bData.surfaceTensionCoeff = sigma;
        bData.EotvosNumber = Eotvos;
        bData.MortonNumber = Morton;

        if( frontMeshInputOption_ == "roughlyCalculated")
        {
           scalar frontMeshLengthScale = 0;
           nps = 0;
           do
           {
               nps = nps + 1;
               scalar bubbleSurfArea = 4.0 * mathematical::pi * pow(radin, 2.0);
               frontMeshLengthScale = pow( 2 * bubbleSurfArea/ ( 2*nps+2*nps*(nps-1)+2*nps*nps*(4-1)), 0.5);
           } while( 0.5*(maxEdge_+minEdge_) > frontMeshLengthScale || frontMeshLengthScale > 0.95*maxEdge_ );
           Info << "       *The front mesh LengthScale is calculated as " << frontMeshLengthScale << endl;
        }
        else if( frontMeshInputOption_ == "finerCalculated" || frontMeshInputOption_ == "manually" )
        {
           scalar frontMeshLengthScale = 0;
           nps = 0;
           do
           {
               nps = nps + 1;
               scalar bubbleSurfArea = 4.0 * mathematical::pi * pow(radin, 2.0);
               frontMeshLengthScale = pow( 2 * bubbleSurfArea/ ( 2*nps+2*nps*(nps-1)+2*nps*nps*(4-1)), 0.5);
           } while( minEdge_ > frontMeshLengthScale || frontMeshLengthScale > (1.0+howFinner_)*minEdge_ );
           Info << "       *The front mesh LengthScale is calculated as " << frontMeshLengthScale << endl;
        }

        scalar dph = mathematical::piByTwo/double(nps);
        label nptot = 2+(4-1)*nps*nps+(nps-1)*nps+nps;
        label netot = 2*nps+2*nps*(nps-1)+2*nps*nps*(4-1);

        Info << "\n       The bubble number " << bI << " is constructed from : " << ' ' << infileName << ' '
             << " with radin/nps/bubblepos/e  as : " << ' ' << radin << ' '
             << nps << ' ' << xc << ' ' << yc << ' '
             << zc  << token::SPACE << e  << token::SPACE << endl;

         Info << "            the outer fluid and bubble densities are : "
              << bData.outerFluidDensity << ' ' << bData.density << endl;

         Info << "            the outer fluid and bubble viscosities are : "
              << bData.outerFluidViscosity << ' ' << bData.viscosity << endl;

         Info << "            the number of point/element for this bubble is : " << ' ' << nptot << ' '
              << netot << ' ' << endl;

       //cccccccccccccccccccccccccccccccccccccccccccccccccccccc
       //c ee is nonzero to create deformed bubbles
            scalar ee = e;

       //cccccccccccccccccccccccccccccccccccccccccccccccccccccc
       //c set north and south pole----CHECK
       scalar rad = radin-ee;

       scalar x = xc;
       scalar y = yc;
       scalar z = zc - rad;
       ptI.currentPoint = point(x, y, z);
       ptI.posInDomain = ptI.currentPoint;
       ptI.currentIndex = 0;
       ptI.keepPoint = true;
       pL.append(ptI);

       x = xc;
       y = yc;
       z = zc + rad;
       ptI.currentPoint = point(x, y, z);
       ptI.posInDomain = ptI.currentPoint;
       ptI.currentIndex = 1;
       ptI.keepPoint = true;
       pL.append(ptI);

       for(int iq = 1; iq <= 4 ; iq++)
       {//1
         for(int i2 = 1; i2 <= nps ; i2++)
         {//2
           for(int i1 = 1; i1 <= nps ; i1++)
           {//3

	       label iip = 2+(iq-1)*nps*nps+(i2-1)*nps+i1;
	       scalar phi = dph*double(i1-i2);
	       label ist = i2-1;
	       if( (i1-i2) < 0)
               {
                 ist=i1-1;
               }
	       scalar theta = mathematical::piByTwo*( double(iq-1) + double(ist)/( double(nps-abs(i1-i2)) + ROOTVSMALL ) );
	       rad=radin+ee*cos(2.0*phi);

	       x = xc + rad*cos(phi)*cos(theta);
	       y = yc + rad*cos(phi)*sin(theta);
	       z = zc - rad*sin(phi);
	       ptI.currentPoint = point(x, y, z);
               ptI.posInDomain = ptI.currentPoint;
	       ptI.currentIndex = iip - 1;
	       ptI.keepPoint = true;
	       pL.append(ptI);

	       label iie=2*i1+2*nps*(i2-1)+2*nps*nps*(iq-1);
	       label ia=iip;
	       label ib=iip+nps;
	       label ic=ib+1;
	       label id=ia+1;

	       label iqq=0;

	       if(i1 == nps)
               {
	           iqq=iq;
	           if(iqq == 4)
                   {
                      iqq=0;
                   }
		   ic=2+iqq*nps*nps+nps+1-i2;
		   id=ic+1;
	       }

	       if(i2 == nps)
               {
	           iqq=iq;
	           if(iqq == 4)
                   {
                      iqq=0;
                   }
		   ib=2+iqq*nps*nps+(nps+1-i1)*nps+1;
		   ic=ib-nps;
	       }

	       if((i1 == nps) && (i2 == 1))
               {
                  id=1;
               }

	       if((i2 == nps) && (i1 == 1))
               {
                  ib=2;
               }

	       elI.keepElement = true;
	       elI.currentIndex = iie-2;
	       elI.pointIndex[0] = ia-1;
	       elI.pointIndex[1] = ic-1;
	       elI.pointIndex[2] = ib-1;
	       eL.append(elI);

	       elI.keepElement = true;
	       elI.currentIndex = iie-1;
	       elI.pointIndex[0] = ia-1;
	       elI.pointIndex[1] = id-1;
	       elI.pointIndex[2] = ic-1;
	       eL.append(elI);

            }//3
          }//2
        }//1

        startPtIndexFromZero(pL,eL);
        findNighboures(pL,eL);
	    pointToElementMapping(pL,eL);
	    ptConnectedPoints(pL,eL);
        printInitialFronts(printAllbubbles[bI], pL, eL);
        initialBubbleVolume_ = calcBubbleVolume(eL);//it's new
        bDataL_.append(bData);
    }//0
    initialAllocatingTheCPTList();
    allCellsInEeachMasket_.setSize(bDataL_.size());
  //allCellsNearEeachFront_.setSize(bDataL_.size());
    setLengthScalesNearTheFront();
}

//--------------------------------------------------------------------------
//-----------------------------findNighboures-------------------------------
//--------------------------------------------------------------------------
template<class CloudType>
inline void Foam::FrontDynamic<CloudType>::findNighboures(DynamicList<pointData>& pL, DynamicList<elementInfo>& eL)
{
    forAll(eL , fI)
    {//1
        elementInfo& elfI = eL[fI];
        elfI.elementIndex[0] = -1;//JJ
        elfI.elementIndex[1] = -1;//JJ
        elfI.elementIndex[2] = -1;//JJ

        label pfI = elfI.pointIndex[0];
        label pfII = elfI.pointIndex[1];
        label pfIII = elfI.pointIndex[2];

        forAll(eL , sI)
        {//2
            elementInfo elsI = eL[sI];

            label psI = elsI.pointIndex[0];
            label psII = elsI.pointIndex[1];
            label psIII = elsI.pointIndex[2];

            if( sI != fI)
            {
                if(  (psII == pfI && psI == pfII) || (psIII == pfI && psII == pfII) || (psI == pfI && psIII == pfII)  )
                {
                    elfI.elementIndex[0] = elsI.currentIndex;
                }
                if(  (psII == pfII && psI == pfIII) || (psIII == pfII && psII == pfIII) || (psI == pfII && psIII == pfIII)  )
                {
                    elfI.elementIndex[1] = elsI.currentIndex;
                }
                if(  (psII == pfIII && psI == pfI) || (psIII == pfIII && psII == pfI) || (psI == pfIII && psIII == pfI)  )
                {
                    elfI.elementIndex[2] = elsI.currentIndex;
                }
            }
        }//2
    }//1
}
//--------------------------------------------------------------------------
//-----------------------------pointToElementMapping---------------------------
//--------------------------------------------------------------------------
// note: func to find connected points of each point in pL.
template<class CloudType>
inline void Foam::FrontDynamic<CloudType>::pointToElementMapping(DynamicList<pointData>& pL, DynamicList<elementInfo>& eL)
{
    forAll(eL, II)
    {//4
        elementInfo& elI = eL[II];

        label p0 = elI.pointIndex[0];
        label p1 = elI.pointIndex[1];
        label p2 = elI.pointIndex[2];

        elI.points[0] = pL[p0].currentPoint;
        elI.points[1] = pL[p1].currentPoint;
        elI.points[2] = pL[p2].currentPoint;

        pL[p0].elOwner = elI.currentIndex;
        pL[p1].elOwner = elI.currentIndex;
        pL[p2].elOwner = elI.currentIndex;
    }//4
}
//--------------------------------------------------------------------------
//-----------------------------ptConnectedPoints---------------------------
//--------------------------------------------------------------------------
// note: func to find connected points of each point in pL.
template<class CloudType>
inline void Foam::FrontDynamic<CloudType>::ptConnectedPoints( DynamicList<pointData>& pL, DynamicList<elementInfo>& eL)
{
    label posInList;
    label f, cNI, cNII, cNIII;
    label curPI, firstPI;

    forAll(pL, I)
    {//1
        pointData& currentPt = pL[I];
        DynamicList<point>& cPoints = currentPt.connectedPoints;
        DynamicList<label>& cPointsIndex = currentPt.connectedPointsIndex;

        cPoints.clear();
        cPointsIndex.clear();
        //cPoints.clearStorage();//JJ
        //cPointsIndex.clearStorage();//JJ

        curPI = currentPt.currentIndex;
        elementInfo elTemp = eL[currentPt.elOwner];

        for( f=0; f<=3; f++)
        {
            if( elTemp.pointIndex[f] == curPI )
            {
                break;
            }
            if( f == 3)
            {
		FatalErrorIn( " ^^^^^^^^^^^^^ we can not find such a point, curPI.") << abort(FatalError);
            }
        }
        findPointsOrder( f, cNI, cNII, cNIII );

        firstPI = elTemp.pointIndex[cNII];
        cPoints.append( elTemp.points[cNII]);
        cPointsIndex.append( elTemp.pointIndex[cNII]);

        label counter = 0;
        do
        {
            elTemp = eL[elTemp.elementIndex[cNIII]];

            for( f=0; f<=3; f++)
            {
                if( elTemp.pointIndex[f] == curPI )
                {
                    break;
                }
                if( f == 3)
                {
		    FatalErrorIn( " ^^^^^^^^^^^^^ we can not find such a point, curPI.") << abort(FatalError);
                }
            }
            findPointsOrder(f, cNI, cNII, cNIII );

            cPoints.append( elTemp.points[cNII]);
            cPointsIndex.append( elTemp.pointIndex[cNII]);

            counter++;
            if( counter > 100)
            {
		FatalErrorIn( " >>>>>>>>>>>>   Sorry, some thing are going to be wrong with You !!!! ")
                << "\n point index is: " << curPI
                << "\n el owner is: " << currentPt.elOwner
                << "\n el owner neighboures: " << eL[currentPt.elOwner].elementIndex[0] << ' '
                                               << eL[currentPt.elOwner].elementIndex[1] << ' '
                                               << eL[currentPt.elOwner].elementIndex[2] << abort(FatalError);
            }
        }while( elTemp.pointIndex[cNIII] != firstPI );
    }//1
}

//--------------------------------------------------------------------------
//-----------------------------printInitialFronts---------------------------
//--------------------------------------------------------------------------
//--------------THIS IS NOT ACCEBTABLE FIX THE WAY IT WRITES TO FILE--------
//- note : func to print the front mesh in .plt format
//- there is an important point that the number
//- of points should be strated from 1  in plt format
template<class CloudType>
inline void Foam::FrontDynamic<CloudType>::printInitialFronts(const fileName& outfileName, DynamicList<pointData> pL, DynamicList<elementInfo> eL)
{
    ofstream outFile(outfileName.c_str(), ios::out);
    ostream& os = outFile;

    // Write header
    os  << " VARIABLES = \"X\" \"Y\" \"Z\",\"NUM\" "
	<< "\n ZONE T=\"Bubble\" "
	<< "\n N=        " << pL.size()
	<< " E=        " << eL.size()
	<< " ZONETYPE=FETriangle"
	<< "\n DATAPACKING=POINT \n";

    // Write vertex coords
    forAll(pL, pI)
    {
        vector curPoint = pL[pI].currentPoint;
	os  << curPoint.x() << token::SPACE
	    << curPoint.y() << token::SPACE
	    << curPoint.z() << token::SPACE << 1 << " \n";
    }

    forAll(eL, eI)
    {
	 os  << eL[eI].pointIndex[0]+1 << token::SPACE
	     << eL[eI].pointIndex[1]+1 << token::SPACE
	     << eL[eI].pointIndex[2]+1 << token::SPACE << " \n";
    }

    outFile.close();
}








//--------------------------------------------------------------------------
//-----------------------------readInitialBubbles---------------------------
//--------------------------------------------------------------------------
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
