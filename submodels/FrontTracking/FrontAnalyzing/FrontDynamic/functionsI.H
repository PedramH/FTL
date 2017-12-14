// note: this cloud cleaning function
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

// note: this func will insert all points
// according to the pointData list entring as an argument.
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

//note: func to calculate the posInDomain of all points and element centres.
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

//note: this function tranfers data from parcel to front.
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

// note: func to arrange the nighboures of every element
// in the eL and the conneceted points for -
// every points in pL.
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

//- note: func to find cNI, cNII, and cNIII, according to point label of curP for current pointData.
template<class CloudType>
inline void Foam::FrontDynamic<CloudType>::findPointsOrder( label f, label& cNI, label& cNII, label& cNIII )
{
    cNI = f;
    cNII = cNI +1;
    if( cNII > 2)
    {
        cNII = cNII -3;
    }
    cNIII = cNII +1;
    if( cNIII > 2)
    {
        cNIII = cNIII - 3;
    }
}


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

//- note: does it start point index from zero? if yes, this returns true.
template<class CloudType>
inline void Foam::FrontDynamic<CloudType>::startPtIndexFromZero(DynamicList<pointData>& pL, DynamicList<elementInfo>& eL)
{
    forAll(pL, pI)
    {    
        if ( pL[pI].currentIndex == pL.size() )
        {
            forAll(pL, pII)
            {
               pL[pII].currentIndex = pL[pII].currentIndex - 1;
               Info << "\n XXXXXXXXXXXXXXXXXXXXXXXX" << endl;
            }
            break;
        }
    }
    forAll(eL, eI)
    {
        label pEI = eL[eI].pointIndex[0];
        label pEII = eL[eI].pointIndex[1];
        label pEIII = eL[eI].pointIndex[2];
	if ( ( pEI == pL.size() ) || ( pEII == pL.size() ) || ( pEIII == pL.size() ) )
        {
               Info << "\n XXXXXXXXXXXXXXXXXXXXXXXX" << endl;
            forAll(eL, eII)
            {
                eL[eII].pointIndex[0] = eL[eII].pointIndex[0] - 1;
                eL[eII].pointIndex[1] = eL[eII].pointIndex[1] - 1;
                eL[eII].pointIndex[2] = eL[eII].pointIndex[2] - 1;
            }
            break;
        }
    }
}    

//- note: func to detect small element. order the lengths, s1 shortest, s3 longest 
template<class CloudType>
inline bool Foam::FrontDynamic<CloudType>::isElementSmall(elementInfo el, scalar minEdge, scalar maxEdge,
                                                          scalar maxAspectRatio, scalar& s1, scalar& s2,
                                                          scalar& s3, label& n1, label& n2, label& n3)
{
    scalar ss = 0;
    label nn = 0;

    if( s1 >= s2)
    {
        ss = s1;
        s1 = s2;
        s2 = ss;
        nn = n1;
        n1 = n2;
        n2 = nn;
    }
     if( s2 >= s3)
    {
        ss = s2;
        s2 = s3;
        s3 = ss;
        nn = n2;
        n2 = n3;
        n3 = nn;
    }
     if( s1 >= s2)
    {
        ss = s1;
        s1 = s2;
        s2 = ss;
        nn = n1;
        n1 = n2;
        n2 = nn;
    }

     if( s1 < minEdge )
     {
        return true;
     }

    return false;
}

//-note: func to detect an element required to be refined. order the lengths, s1 shortest, s3 longest 
template<class CloudType>
inline bool Foam::FrontDynamic<CloudType>::isElementRefined(elementInfo el, scalar minEdge, scalar maxEdge, scalar maxAspectRatio,
                                                            scalar& s1, scalar& s2, scalar& s3, label& n1, label& n2, label& n3)
{
    vector e01 = ( el.points[1] - el.points[0]);
    vector e02 = ( el.points[2] - el.points[0]);

    scalar ss;//JJ
    label nn;//JJ

    if( s1 >= s2)
    {
        ss = s1;
        s1 = s2;
        s2 = ss;
        nn = n1;
        n1 = n2;
        n2 = nn;
    }
    if( s2 >= s3)
    {
        ss = s2;
        s2 = s3;
        s3 = ss;
        nn = n2;
        n2 = n3;
        n3 = nn;
    }
    if( s1 >= s2)
    {
        ss = s1;
        s1 = s2;
        s2 = ss;
        nn = n1;
        n1 = n2;
        n2 = nn;
    }

    scalar surfaceArea = mag (0.5 * (e01 ^ e02));
    //scalar elAspectRatio = 0.25 * sqrt (3.0) * sqr( ( s1 + s2 + s3) / 3.) * surfaceArea;
    scalar elAspectRatio = sqr(s1 + s2 + s3)/( sqrt (3.0) * 12.0 * surfaceArea);


    if( s3 > maxEdge )
    {
        return true;
    }
    
    if( s1 > minEdge && elAspectRatio > maxAspectRatio )
    {
        //Info << " aspect rario " << elAspectRatio << endl;
        return true;
    }
    
    return false;
}


//- note: func to coarse the surface grid
template<class CloudType>
inline void Foam::FrontDynamic<CloudType>::frontCoarsening(DynamicList<pointData>& pL, DynamicList<elementInfo>& eL,
                                                           scalar minEdge, scalar maxEdge, scalar maxAspectRatio)
{//1
    label pI1, pI2, pI3, pI4;
    label c4I;
    label c4II;
    label c4III;
    label c5I;
    label c5II;
    label c5III;
    label c6I;
    label c6II;
    label c6III;
    label c7I;
    label c7II;
    label c7III;
    label c8I;
    label c8II;
    label c8III;
    label c9I, c9II, c9III;
    label c10I, c10II, c10III;
    label c11I, c11II, c11III;
    label I, II, III;
    label f, cNI, cNII, cNIII;

    label numberOfELDeleted = 0;

    label initialElListSize = eL.size();

    forAll(eL, elementI)
    //for(int elementI=0; elementI < initialElListSize; elementI++)   
    //for(int elementI=0; elementI < initial_; elementI++)//JJ   
    {//2

        elementInfo& currentElemet = eL[elementI];

        // order the lengths, s1 shortest, s3 longest 
        scalar s1 = mag( currentElemet.points[1] - currentElemet.points[0]);
        scalar s2 = mag( currentElemet.points[2] - currentElemet.points[1]);
        scalar s3 = mag( currentElemet.points[0] - currentElemet.points[2]);
        label n1 = 0;
        label n2 = 1;
        label n3 = 2;
        bool smallElement = isElementSmall(currentElemet, minEdge, maxEdge, maxAspectRatio, s1, s2, s3, n1, n2, n3);

        //Info<< "Hey, el num: " << ' ' << elementI << endl;
        //Info<< "s1/s2/s3 " << ' ' << s1 << ' ' << s2 << ' ' << s3 << endl;
        //Info<< "AAHey--->smallOk/keppElOk/n1 " << ' ' << smallElement << ' ' << currentElemet.keepElement << ' ' << n1 << endl;

        if ( smallElement && currentElemet.keepElement)
        {//3
            //Info<< "The element " << ' ' << currentElemet.currentIndex << ' ' << " will be coarsed." << endl;

            numberOfELDeleted = numberOfELDeleted + 2;

            label shortEdge = n1;

            // I, II, III are the point indexes of current element.
            // n1 is the index of shortest edge and p1.
            findPointsOrder( shortEdge, I, II, III);
            point p1 = currentElemet.points[I];//JJ
            point p2 = currentElemet.points[II];//JJ
            point p3 = currentElemet.points[III];//JJ
            pI1 = currentElemet.pointIndex[I];
            pI2 = currentElemet.pointIndex[II];
            pI3 = currentElemet.pointIndex[III];

            //Info << "all Nbs---> " << currentElemet.elementIndex[I] << ' ' << currentElemet.elementIndex[II] << ' ' 
            //     << currentElemet.elementIndex[III] << endl;

            // find point 5, third point of neighbElementII
            //label posInList = quickSelect(eL, currentElemet.elementIndex[II], 0, eL.size()-1);
            //elementInfo& neighbElementII = eL[posInList];
            elementInfo& neighbElementII = eL[currentElemet.elementIndex[II]];
            for( f=0; f<=3; f++)
            {
                if( neighbElementII.pointIndex[f] == pI2)
                {
                    //Info<< "f5 " << ' ' << f<< endl;
                    break;
                }
                if( f == 3)
                {
		    FatalErrorIn
		    (
			    " ^^^^^^^^^^^^^ we can not find such a point."
		    ) << abort(FatalError);
                }
            }            
            // c5I, c5II, c5III are the point indexes of neighbElementII.
            findPointsOrder( f, cNI, cNII, cNIII);
            point p5 = neighbElementII.points[cNII];
            c5I = cNII;
	    c5II = cNIII;
	    c5III = cNI;

            // find point 6, third point of neighbElementIII
            //posInList = quickSelect(eL, currentElemet.elementIndex[III], 0, eL.size()-1);
            //elementInfo& neighbElementIII = eL[posInList];
            elementInfo& neighbElementIII = eL[currentElemet.elementIndex[III]];
            for( f=0; f<=3; f++)
            {
                if( neighbElementIII.pointIndex[f] == pI3)
                {
                    //Info<< "f6 " << ' ' << f<< endl;
                    break;
                }
                if( f == 3)
                {
		    FatalErrorIn
		    (
			    " ^^^^^^^^^^^^^ we can not find such a point."
		    ) << abort(FatalError);
                }
            }
                    //Info<< "f6 " << ' ' << f<< endl;
            findPointsOrder( f, cNI, cNII, cNIII);
            point p6 = neighbElementIII.points[cNII];
	    c6I = cNII;
	    c6II = cNIII;
	    c6III = cNI;

            // find point 4, third point of neighbElementIII
            //posInList = quickSelect(eL, currentElemet.elementIndex[I], 0, eL.size()-1);
            //elementInfo& neighbElementI = eL[posInList];
            elementInfo& neighbElementI = eL[currentElemet.elementIndex[I]];
            for( f=0; f<=3; f++)
            {
                if( neighbElementI.pointIndex[f] == pI1)
                {
                    //Info<< "f4 " << ' ' << f<< endl;
                    break;
                }
                if( f == 3)
                {
		    FatalErrorIn
		    (
			    " ^^^^^^^^^^^^^ we can not find such a point."
		    ) << abort(FatalError);
                }
            }
            findPointsOrder( f, cNI, cNII, cNIII);
            point p4 = neighbElementI.points[cNII];
            pI4 = neighbElementI.pointIndex[cNII];
	    c4I = cNII;
	    c4II = cNIII;
	    c4III = cNI;


            // find point 7, third point of neighbElementILeft
            //posInList = quickSelect(eL, neighbElementI.elementIndex[c4III], 0, eL.size()-1);
            //elementInfo& neighbElementILeft = eL[posInList];
            elementInfo& neighbElementILeft = eL[neighbElementI.elementIndex[c4III]];
            for( f=0; f<=3; f++)
            {
                if( neighbElementILeft.pointIndex[f] == pI1)
                {
                    //Info<< "f7 " << ' ' << f<< endl;
                    break;
                }
                if( f == 3)
                {
		    FatalErrorIn
		    (
			    " ^^^^^^^^^^^^^ we can not find such a point."
		    ) << abort(FatalError);
                }
            }
            findPointsOrder( f, cNI, cNII, cNIII);
            point p7 = neighbElementILeft.points[cNII];
	    c7I = cNII;
	    c7II = cNIII;
	    c7III = cNI;

            // find point 8, third point of neighbElementIRight
            //posInList = quickSelect(eL, neighbElementI.elementIndex[c4I], 0, eL.size()-1);
            //elementInfo& neighbElementIRight = eL[posInList];
            elementInfo& neighbElementIRight = eL[neighbElementI.elementIndex[c4I]];
            for( f=0; f<=3; f++)
            {
                if( neighbElementIRight.pointIndex[f] == pI4)
                {
                    //Info<< "f8 " << ' ' << f<< endl;
                    break;
                }
                if( f == 3)
                {
		    FatalErrorIn
		    (
			    " ^^^^^^^^^^^^^ we can not find such a point."
		    ) << abort(FatalError);
                }
            }
            findPointsOrder( f, cNI, cNII, cNIII);
            point p8 = neighbElementIRight.points[cNII];
	    c8I = cNII;
	    c8II = cNIII;
	    c8III = cNI;

            if( neighbElementIRight.currentIndex == neighbElementII.currentIndex )
            {

                Info << " UUUUUUUUUUUUUUUUUUUUUUUUU---the right triple is occured." << currentElemet.currentIndex << endl;

                // allocation of the new point: 
                point newPoint = 0.5 * ( p1 + p2) + 0.125 * ( p3 + p4) - 0.0625 * ( p5 + p6 + p7 + p8);

                // position of p1 now is updated.
                //posInList = quickSelect(pL, pI1, 0, pL.size()-1);
                //pL[posInList].currentPoint = newPoint;
                pL[pI1].currentPoint = newPoint;

                // transform only all positions of points of the element sharing the lable of p1 to newPoint 
                //cNI = c6I;
                //cNII = c6II;
                //cNIII = c6III;
                cNIII = c6I; 
                cNI = c6II;                           
                elementInfo elTemp1 = neighbElementIII;
                neighbElementIII.points[cNI] = newPoint;
                do 
                {
                    //label posInList = quickSelect(eL, elTemp1.elementIndex[cNI], 0, eL.size()-1);
                    //elementInfo& neighbElTemp1 = eL[posInList];
                    elementInfo& neighbElTemp1 = eL[elTemp1.elementIndex[cNIII]];
                    for( f=0; f<=3; f++)
                    {
                       if( neighbElTemp1.pointIndex[f] == pI1 )
                       {
                           //Info<< "F I can break " << ' ' << f << endl;
                           break;
                       }
                       if( f == 3)
                       {
		           FatalErrorIn
		           (
			       " ^^^^^^^^^^^^^ we can not find such a point, pI1."
		           ) << abort(FatalError);
                       }
                    }
                    findPointsOrder( f, cNI, cNII, cNIII);                        
                    neighbElTemp1.points[cNI] = newPoint;
                    elTemp1 = neighbElTemp1;
                } while( elTemp1.currentIndex != neighbElementILeft.currentIndex );

                // transform all positions and indexes of points of the element sharing the p2 to newPoint
                // then the p2 must be eliminated.
                neighbElementII.points[c5III] = newPoint;
                neighbElementII.pointIndex[c5III] = pI1;

                // update the neighboures of neighbElementIII and neighbElementII 
                neighbElementIII.elementIndex[c6II] = neighbElementII.currentIndex;
                neighbElementII.elementIndex[c5II] = neighbElementIII.currentIndex;

                // update the neighboures of neighbElementILeft and neighbElementIRight 
                neighbElementILeft.elementIndex[c7II] = neighbElementIRight.currentIndex;
                neighbElementIRight.elementIndex[c8II] = neighbElementILeft.currentIndex;                 

                // p2 now is eliminated.
                //posInList = quickSelect(pL, pI2, 0, pL.size()-1);
                //pL[posInList].keepPoint = false;
                pL[pI2].keepPoint = false;

                // eleminating the elements: currentElemet and neighbElementI
                currentElemet.keepElement = false;
                neighbElementI.keepElement = false;

            }
            else if( neighbElementILeft.currentIndex == neighbElementIII.currentIndex  )
            {

                Info << " UUUUUUUUUUUUUUUUUUUUUUUUU---the left triple is occured." << currentElemet.currentIndex << endl;

                // allocation of the new point: 
                point newPoint = 0.5 * ( p1 + p2) + 0.125 * ( p3 + p4) - 0.0625 * ( p5 + p6 + p7 + p8);

                // position of p1 now is updated.
                //posInList = quickSelect(pL, pI1, 0, pL.size()-1);
                //pL[posInList].currentPoint = newPoint;
                pL[pI1].currentPoint = newPoint;

                // transform only all positions of points of the element sharing the lable of p1 to newPoint 
                neighbElementIII.points[c6II] = newPoint;

                // transform all positions and indexes of points of the element sharing the p2 to newPoint
                // then the p2 must be eliminated.
                //cNI = c8I;
                //cNII = c8II;
                //cNIII = c8III;             
                cNIII = c8I;
                cNI = c8II;                          
                elementInfo elTemp2 = neighbElementIRight;
                neighbElementIRight.points[cNI] = newPoint;
                neighbElementIRight.pointIndex[cNI] = pI1;
                do 
                {
                    //label posInList = quickSelect(eL, elTemp2.elementIndex[cNI], 0, eL.size()-1);
                    //elementInfo& neighbElTemp2 = eL[posInList];
                    elementInfo& neighbElTemp2 = eL[elTemp2.elementIndex[cNIII]];
                    for( f=0; f<=3; f++)
                    {
                        if( neighbElTemp2.pointIndex[f] == pI2 )
                        {
                            //Info<< "S I can break " << ' ' << f << endl;
                            break;
                        }
                        if( f == 3)
                        {
		            FatalErrorIn
		            (
			        " ^^^^^^^^^^^^^ we can not find such a point, pI2."
		            ) << abort(FatalError);
                        }
                    }
                    findPointsOrder( f, cNI, cNII, cNIII);
                    neighbElTemp2.points[cNI] = newPoint;
                    neighbElTemp2.pointIndex[cNI] = pI1;
                    elTemp2 = neighbElTemp2;
                } while( elTemp2.currentIndex != neighbElementII.currentIndex );
                                      
                // update the neighboures of neighbElementIII and neighbElementII 
                neighbElementIII.elementIndex[c6II] = neighbElementII.currentIndex;
                neighbElementII.elementIndex[c5II] = neighbElementIII.currentIndex;

                // update the neighboures of neighbElementILeft and neighbElementIRight 
                neighbElementILeft.elementIndex[c7II] = neighbElementIRight.currentIndex;
                neighbElementIRight.elementIndex[c8II] = neighbElementILeft.currentIndex;                 

                // p2 now is eliminated.
                //posInList = quickSelect(pL, pI2, 0, pL.size()-1);
                //pL[posInList].keepPoint = false;
                pL[pI2].keepPoint = false;

                // eleminating the elements: currentElemet and neighbElementI
                currentElemet.keepElement = false;
                neighbElementI.keepElement = false;
            }

            else
            {
                // allocation of the new point: 
                point newPoint = 0.5 * ( p1 + p2) + 0.125 * ( p3 + p4) - 0.0625 * ( p5 + p6 + p7 + p8);

                // position of p1 now is updated.
                //posInList = quickSelect(pL, pI1, 0, pL.size()-1);
                //pL[posInList].currentPoint = newPoint;
                pL[pI1].currentPoint = newPoint;

                // transform only all positions of points of the element sharing the lable of p1 to newPoint 
                //cNI = c6I;
                //cNII = c6II;
                //cNIII = c6III;
                cNIII = c6I; 
                cNI = c6II;                           
                elementInfo elTemp1 = neighbElementIII;
                neighbElementIII.points[cNI] = newPoint;
                do 
                {
                    //label posInList = quickSelect(eL, elTemp1.elementIndex[cNI], 0, eL.size()-1);
                    //elementInfo& neighbElTemp1 = eL[posInList];
                    elementInfo& neighbElTemp1 = eL[elTemp1.elementIndex[cNIII]];
                    for( f=0; f<=3; f++)
                    {
                       if( neighbElTemp1.pointIndex[f] == pI1 )
                       {
                           //Info<< "F I can break " << ' ' << f << endl;
                           break;
                       }
                       if( f == 3)
                       {
		           FatalErrorIn
		           (
			       " ^^^^^^^^^^^^^ we can not find such a point, pI1."
		           ) << abort(FatalError);
                       }
                    }
                    findPointsOrder( f, cNI, cNII, cNIII);                        
                    neighbElTemp1.points[cNI] = newPoint;
                    elTemp1 = neighbElTemp1;
                         //Info<< "F elTemp1.currentIndex/neighbElementILeft.currentIndex " << elTemp1.currentIndex << ' ' << neighbElementILeft.currentIndex << endl;
                } while( elTemp1.currentIndex != neighbElementILeft.currentIndex );

                // transform all positions and indexes of points of the element sharing the p2 to newPoint
                // then the p2 must be eliminated.
                //cNI = c8I;
                //cNII = c8II;
                //cNIII = c8III;             
                cNIII = c8I;
                cNI = c8II;                          
                elementInfo elTemp2 = neighbElementIRight;
                neighbElementIRight.points[cNI] = newPoint;
                neighbElementIRight.pointIndex[cNI] = pI1;

                         //Info << "start/neighbour/currentElemet.currentIndex  " << elTemp2.currentIndex << ' '
                         //     << elTemp2.elementIndex[cNIII] << ' ' << currentElemet.currentIndex << endl;
                         //Info << "el/II/III/I/left/right  " << currentElemet.currentIndex << ' '
                         //     << currentElemet.elementIndex[II] << ' ' << currentElemet.elementIndex[III] << ' ' << currentElemet.elementIndex[I] << ' ' 
                         //     << neighbElementI.elementIndex[c4III] << ' ' << neighbElementI.elementIndex[c4I] << endl;

                         //Info << "p1/p2/p3/p4/p5/p6/p7/p8  " << pI1 << ' ' << pI2 << ' ' << pI3 << ' ' << pI4 << ' ' 
                         //     << neighbElementII.pointIndex[c5I] << ' ' << neighbElementIII.pointIndex[c6I] << ' ' 
                         //     << neighbElementILeft.pointIndex[c7I] << ' ' << neighbElementIRight.pointIndex[c8I] << endl;
                do 
                {
                    //label posInList = quickSelect(eL, elTemp2.elementIndex[cNI], 0, eL.size()-1);
                    //elementInfo& neighbElTemp2 = eL[posInList];
                    elementInfo& neighbElTemp2 = eL[elTemp2.elementIndex[cNIII]];
                    for( f=0; f<=3; f++)
                    {
                        if( neighbElTemp2.pointIndex[f] == pI2 )
                        {
                            //Info<< "S I can break " << ' ' << f << endl;
                            break;
                        }
                        if( f == 3)
                        {
		            FatalErrorIn
		            (
			        " ^^^^^^^^^^^^^ we can not find such a point, pI2."
		            ) << abort(FatalError);
                        }
                    }
                    findPointsOrder( f, cNI, cNII, cNIII);
                    neighbElTemp2.points[cNI] = newPoint;
                    neighbElTemp2.pointIndex[cNI] = pI1;
                    elTemp2 = neighbElTemp2;
                         //Info<< "S elTemp2.currentIndex/neighbElementII.currentIndex " << elTemp2.currentIndex << ' ' << neighbElementII.currentIndex << endl;
                } while( elTemp2.currentIndex != neighbElementII.currentIndex );
                                      
                // update the neighboures of neighbElementIII and neighbElementII 
                neighbElementIII.elementIndex[c6II] = neighbElementII.currentIndex;
                neighbElementII.elementIndex[c5II] = neighbElementIII.currentIndex;

                // update the neighboures of neighbElementILeft and neighbElementIRight 
                neighbElementILeft.elementIndex[c7II] = neighbElementIRight.currentIndex;
                neighbElementIRight.elementIndex[c8II] = neighbElementILeft.currentIndex;                 

                // p2 now is eliminated.
                //posInList = quickSelect(pL, pI2, 0, pL.size()-1);
                //pL[posInList].keepPoint = false;
                pL[pI2].keepPoint = false;

                // eleminating the elements: currentElemet and neighbElementI
                currentElemet.keepElement = false;
                neighbElementI.keepElement = false;

                // updatong the all elOwner
                //posInList = quickSelect(pL, pI1, 0, pL.size()-1);
                //pL[posInList].elOwner = neighbElementILeft.currentIndex;
                //pL[pI3].elOwner = neighbElementIII.currentIndex;
                //pL[pI1].elOwner = neighbElementIII.currentIndex;
                //pL[pI4].elOwner = neighbElementIRight.currentIndex;
            }
        }//3
                 
    }//2
   Info << "                 After coarsening " << numberOfELDeleted << " element/elemens were deleted." << endl;
}//1

/*

//- note: func to coarse the surface grid
template<class CloudType>
inline void Foam::FrontDynamic<CloudType>::frontCoarsening(DynamicList<pointData>& pL, DynamicList<elementInfo>& eL,
                                                           scalar minEdge, scalar maxEdge, scalar maxAspectRatio)
{//1
    label pI1, pI2, pI3, pI4;
    label c4I;
    label c4II;
    label c4III;
    label c5I;
    label c5II;
    label c5III;
    label c6I;
    label c6II;
    label c6III;
    label c7I;
    label c7II;
    label c7III;
    label c8I;
    label c8II;
    label c8III;
    label c9I, c9II, c9III;
    label c10I, c10II, c10III;
    label c11I, c11II, c11III;
    label I, II, III;
    label f, cNI, cNII, cNIII;

    label numberOfELDeleted = 0;

    label initialElListSize = eL.size();

    forAll(eL, elementI)
    //for(int elementI=0; elementI < initialElListSize; elementI++)   
    //for(int elementI=0; elementI < initial_; elementI++)//JJ   
    {//2

        elementInfo& currentElemet = eL[elementI];

        // order the lengths, s1 shortest, s3 longest 
        scalar s1 = mag( currentElemet.points[1] - currentElemet.points[0]);
        scalar s2 = mag( currentElemet.points[2] - currentElemet.points[1]);
        scalar s3 = mag( currentElemet.points[0] - currentElemet.points[2]);
        label n1 = 0;
        label n2 = 1;
        label n3 = 2;
        bool smallElement = isElementSmall(currentElemet, minEdge, maxEdge, maxAspectRatio, s1, s2, s3, n1, n2, n3);

        //Info<< "Hey, el num: " << ' ' << elementI << endl;
        //Info<< "s1/s2/s3 " << ' ' << s1 << ' ' << s2 << ' ' << s3 << endl;
        //Info<< "AAHey--->smallOk/keppElOk/n1 " << ' ' << smallElement << ' ' << currentElemet.keepElement << ' ' << n1 << endl;

        if ( smallElement && currentElemet.keepElement)
        {//3
            //Info<< "The element " << ' ' << currentElemet.currentIndex << ' ' << " will be coarsed." << endl;

            numberOfELDeleted = numberOfELDeleted + 2;

            label shortEdge = n1;

            // I, II, III are the point indexes of current element.
            // n1 is the index of shortest edge and p1.
            findPointsOrder( shortEdge, I, II, III);
            point p1 = currentElemet.points[I];//JJ
            point p2 = currentElemet.points[II];//JJ
            point p3 = currentElemet.points[III];//JJ
            pI1 = currentElemet.pointIndex[I];
            pI2 = currentElemet.pointIndex[II];
            pI3 = currentElemet.pointIndex[III];

            //Info << "all Nbs---> " << currentElemet.elementIndex[I] << ' ' << currentElemet.elementIndex[II] << ' ' 
            //     << currentElemet.elementIndex[III] << endl;

            // find point 5, third point of neighbElementII
            //label posInList = quickSelect(eL, currentElemet.elementIndex[II], 0, eL.size()-1);
            //elementInfo& neighbElementII = eL[posInList];
            elementInfo& neighbElementII = eL[currentElemet.elementIndex[II]];
            for( f=0; f<=3; f++)
            {
                if( neighbElementII.pointIndex[f] == pI2)
                {
                    //Info<< "f5 " << ' ' << f<< endl;
                    break;
                }
                if( f == 3)
                {
		    FatalErrorIn
		    (
			    " ^^^^^^^^^^^^^ we can not find such a point."
		    ) << abort(FatalError);
                }
            }            
            // c5I, c5II, c5III are the point indexes of neighbElementII.
            findPointsOrder( f, cNI, cNII, cNIII);
            point p5 = neighbElementII.points[cNII];
            c5I = cNII;
	    c5II = cNIII;
	    c5III = cNI;

            // find point 6, third point of neighbElementIII
            //posInList = quickSelect(eL, currentElemet.elementIndex[III], 0, eL.size()-1);
            //elementInfo& neighbElementIII = eL[posInList];
            elementInfo& neighbElementIII = eL[currentElemet.elementIndex[III]];
            for( f=0; f<=3; f++)
            {
                if( neighbElementIII.pointIndex[f] == pI3)
                {
                    //Info<< "f6 " << ' ' << f<< endl;
                    break;
                }
                if( f == 3)
                {
		    FatalErrorIn
		    (
			    " ^^^^^^^^^^^^^ we can not find such a point."
		    ) << abort(FatalError);
                }
            }
                    //Info<< "f6 " << ' ' << f<< endl;
            findPointsOrder( f, cNI, cNII, cNIII);
            point p6 = neighbElementIII.points[cNII];
	    c6I = cNII;
	    c6II = cNIII;
	    c6III = cNI;

            // find point 4, third point of neighbElementIII
            //posInList = quickSelect(eL, currentElemet.elementIndex[I], 0, eL.size()-1);
            //elementInfo& neighbElementI = eL[posInList];
            elementInfo& neighbElementI = eL[currentElemet.elementIndex[I]];
            for( f=0; f<=3; f++)
            {
                if( neighbElementI.pointIndex[f] == pI1)
                {
                    //Info<< "f4 " << ' ' << f<< endl;
                    break;
                }
                if( f == 3)
                {
		    FatalErrorIn
		    (
			    " ^^^^^^^^^^^^^ we can not find such a point."
		    ) << abort(FatalError);
                }
            }
            findPointsOrder( f, cNI, cNII, cNIII);
            point p4 = neighbElementI.points[cNII];
            pI4 = neighbElementI.pointIndex[cNII];
	    c4I = cNII;
	    c4II = cNIII;
	    c4III = cNI;


            // find point 7, third point of neighbElementILeft
            //posInList = quickSelect(eL, neighbElementI.elementIndex[c4III], 0, eL.size()-1);
            //elementInfo& neighbElementILeft = eL[posInList];
            elementInfo& neighbElementILeft = eL[neighbElementI.elementIndex[c4III]];
            for( f=0; f<=3; f++)
            {
                if( neighbElementILeft.pointIndex[f] == pI1)
                {
                    //Info<< "f7 " << ' ' << f<< endl;
                    break;
                }
                if( f == 3)
                {
		    FatalErrorIn
		    (
			    " ^^^^^^^^^^^^^ we can not find such a point."
		    ) << abort(FatalError);
                }
            }
            findPointsOrder( f, cNI, cNII, cNIII);
            point p7 = neighbElementILeft.points[cNII];
	    c7I = cNII;
	    c7II = cNIII;
	    c7III = cNI;

            // find point 8, third point of neighbElementIRight
            //posInList = quickSelect(eL, neighbElementI.elementIndex[c4I], 0, eL.size()-1);
            //elementInfo& neighbElementIRight = eL[posInList];
            elementInfo& neighbElementIRight = eL[neighbElementI.elementIndex[c4I]];
            for( f=0; f<=3; f++)
            {
                if( neighbElementIRight.pointIndex[f] == pI4)
                {
                    //Info<< "f8 " << ' ' << f<< endl;
                    break;
                }
                if( f == 3)
                {
		    FatalErrorIn
		    (
			    " ^^^^^^^^^^^^^ we can not find such a point."
		    ) << abort(FatalError);
                }
            }
            findPointsOrder( f, cNI, cNII, cNIII);
            point p8 = neighbElementIRight.points[cNII];
	    c8I = cNII;
	    c8II = cNIII;
	    c8III = cNI;

            if( neighbElementIRight.currentIndex == neighbElementII.currentIndex )
            {
                Info << " UUUUUUUUUUUUUUUUUUUUUUUUU---the right triple is occured." << currentElemet.currentIndex << endl;
		//        FatalErrorIn
		//        (
		//	    " UUUUUUUUUUUUUUUUUUUUUUUUU---the right triple is occured."
		//        ) << abort(FatalError);
                currentElemet.points[II] = p4;
                currentElemet.pointIndex[II] = pI4;

                // finding neighbElementIIRight
                elementInfo& neighbElementIIRight = eL[neighbElementII.elementIndex[c5I]];
                for( f=0; f<=3; f++)
                {
                    if( neighbElementIIRight.pointIndex[f] == pI4)
                    {
                        break;
                    }
                    if( f == 3)
                    {
		        FatalErrorIn
		        (
			    " ^^^^^^^^^^^^^ we can not find pI4 for neighbElementIIRight."
		        ) << abort(FatalError);
                    }
                }
                findPointsOrder( f, cNI, cNII, cNIII);
	        c9I = cNII;
	        c9II = cNIII;
	        c9III = cNI;
                neighbElementIIRight.elementIndex[c9II] = currentElemet.currentIndex;
                currentElemet.elementIndex[II] = neighbElementIIRight.currentIndex;

                // finding neighbElementIBottom
                elementInfo& neighbElementIBottom = eL[neighbElementI.elementIndex[c4III]];
                for( f=0; f<=3; f++)
                {
                    if( neighbElementIBottom.pointIndex[f] == pI1)
                    {
                        break;
                    }
                    if( f == 3)
                    {
		        FatalErrorIn
		        (
			    " ^^^^^^^^^^^^^ we can not find pI1 for neighbElementIBottom."
		        ) << abort(FatalError);
                    }
                }
                findPointsOrder( f, cNI, cNII, cNIII);
	        c10I = cNII;
	        c10II = cNIII;
	        c10III = cNI;
                neighbElementIBottom.elementIndex[c10II] = currentElemet.currentIndex;
                currentElemet.elementIndex[I] = neighbElementIBottom.currentIndex;

                pL[pI2].keepPoint = false;
                neighbElementI.keepElement = false;
                neighbElementII.keepElement = false;
            }
            else if( neighbElementILeft.currentIndex == neighbElementIII.currentIndex  )
            {
                Info << " UUUUUUUUUUUUUUUUUUUUUUUUU---the left triple is occured." << currentElemet.currentIndex << endl;
		//        FatalErrorIn
		//        (
		//	    " UUUUUUUUUUUUUUUUUUUUUUUUU---the left triple is occured."
		//        ) << abort(FatalError);
                currentElemet.points[I] = p4;
                currentElemet.pointIndex[I] = pI4;

                // finding neighbElementIIILeft
                elementInfo& neighbElementIIILeft = eL[neighbElementIII.elementIndex[c6III]];
                for( f=0; f<=3; f++)
                {
                    if( neighbElementIIILeft.pointIndex[f] == pI3)
                    {
                        break;
                    }
                    if( f == 3)
                    {
		        FatalErrorIn
		        (
			    " ^^^^^^^^^^^^^ we can not find pI3 neighbElementIIILeft."
		        ) << abort(FatalError);
                    }
                }
                findPointsOrder( f, cNI, cNII, cNIII);
	        c11I = cNII;
	        c11II = cNIII;
	        c11III = cNI;
                neighbElementIIILeft.elementIndex[c11II] = currentElemet.currentIndex;
                currentElemet.elementIndex[III] = neighbElementIIILeft.currentIndex;

                // finding neighbElementIBottom
                elementInfo& neighbElementIBottom = eL[neighbElementI.elementIndex[c4I]];
                for( f=0; f<=3; f++)
                {
                    if( neighbElementIBottom.pointIndex[f] == pI4)
                    {
                        break;
                    }
                    if( f == 3)
                    {
		        FatalErrorIn
		        (
			    " ^^^^^^^^^^^^^ we can not find pI4 for neighbElementIBottom."
		        ) << abort(FatalError);
                    }
                }
                findPointsOrder( f, cNI, cNII, cNIII);
	        c10I = cNII;
	        c10II = cNIII;
	        c10III = cNI;
                neighbElementIBottom.elementIndex[c10II] = currentElemet.currentIndex;
                currentElemet.elementIndex[I] = neighbElementIBottom.currentIndex;

                pL[pI1].keepPoint = false;
                neighbElementI.keepElement = false;
                neighbElementIII.keepElement = false;
            }
            else
            {
                // allocation of the new point: 
                point newPoint = 0.5 * ( p1 + p2) + 0.125 * ( p3 + p4) - 0.0625 * ( p5 + p6 + p7 + p8);

                // position of p1 now is updated.
                //posInList = quickSelect(pL, pI1, 0, pL.size()-1);
                //pL[posInList].currentPoint = newPoint;
                pL[pI1].currentPoint = newPoint;

                // transform only all positions of points of the element sharing the lable of p1 to newPoint 
                //cNI = c6I;
                //cNII = c6II;
                //cNIII = c6III;
                cNIII = c6I; 
                cNI = c6II;                           
                elementInfo elTemp1 = neighbElementIII;
                neighbElementIII.points[cNI] = newPoint;
                do 
                {
                    //label posInList = quickSelect(eL, elTemp1.elementIndex[cNI], 0, eL.size()-1);
                    //elementInfo& neighbElTemp1 = eL[posInList];
                    elementInfo& neighbElTemp1 = eL[elTemp1.elementIndex[cNIII]];
                    for( f=0; f<=3; f++)
                    {
                       if( neighbElTemp1.pointIndex[f] == pI1 )
                       {
                           //Info<< "F I can break " << ' ' << f << endl;
                           break;
                       }
                       if( f == 3)
                       {
		           FatalErrorIn
		           (
			       " ^^^^^^^^^^^^^ we can not find such a point, pI1."
		           ) << abort(FatalError);
                       }
                    }
                    findPointsOrder( f, cNI, cNII, cNIII);                        
                    neighbElTemp1.points[cNI] = newPoint;
                    elTemp1 = neighbElTemp1;
                         //Info<< "F elTemp1.currentIndex/neighbElementILeft.currentIndex " << elTemp1.currentIndex << ' ' << neighbElementILeft.currentIndex << endl;
                } while( elTemp1.currentIndex != neighbElementILeft.currentIndex );

                // transform all positions and indexes of points of the element  sharing the p2 to newPoint
                // then the p2 must be eliminated.
                //cNI = c8I;
                //cNII = c8II;
                //cNIII = c8III;             
                cNIII = c8I;
                cNI = c8II;                          
                elementInfo elTemp2 = neighbElementIRight;
                neighbElementIRight.points[cNI] = newPoint;
                neighbElementIRight.pointIndex[cNI] = pI1;

                         //Info << "start/neighbour/currentElemet.currentIndex  " << elTemp2.currentIndex << ' '
                         //     << elTemp2.elementIndex[cNIII] << ' ' << currentElemet.currentIndex << endl;
                         //Info << "el/II/III/I/left/right  " << currentElemet.currentIndex << ' '
                         //     << currentElemet.elementIndex[II] << ' ' << currentElemet.elementIndex[III] << ' ' << currentElemet.elementIndex[I] << ' ' 
                         //     << neighbElementI.elementIndex[c4III] << ' ' << neighbElementI.elementIndex[c4I] << endl;

                         //Info << "p1/p2/p3/p4/p5/p6/p7/p8  " << pI1 << ' ' << pI2 << ' ' << pI3 << ' ' << pI4 << ' ' 
                         //     << neighbElementII.pointIndex[c5I] << ' ' << neighbElementIII.pointIndex[c6I] << ' ' 
                         //     << neighbElementILeft.pointIndex[c7I] << ' ' << neighbElementIRight.pointIndex[c8I] << endl;
                do 
                {
                    //label posInList = quickSelect(eL, elTemp2.elementIndex[cNI], 0, eL.size()-1);
                    //elementInfo& neighbElTemp2 = eL[posInList];
                    elementInfo& neighbElTemp2 = eL[elTemp2.elementIndex[cNIII]];
                    for( f=0; f<=3; f++)
                    {
                        if( neighbElTemp2.pointIndex[f] == pI2 )
                        {
                            //Info<< "S I can break " << ' ' << f << endl;
                            break;
                        }
                        if( f == 3)
                        {
		            FatalErrorIn
		            (
			        " ^^^^^^^^^^^^^ we can not find such a point, pI2."
		            ) << abort(FatalError);
                        }
                    }
                    findPointsOrder( f, cNI, cNII, cNIII);
                    neighbElTemp2.points[cNI] = newPoint;
                    neighbElTemp2.pointIndex[cNI] = pI1;
                    elTemp2 = neighbElTemp2;
                         //Info<< "S elTemp2.currentIndex/neighbElementII.currentIndex " << elTemp2.currentIndex << ' ' << neighbElementII.currentIndex << endl;
                } while( elTemp2.currentIndex != neighbElementII.currentIndex );
                                      
                // update the neighboures of neighbElementIII and neighbElementII 
                neighbElementIII.elementIndex[c6II] = neighbElementII.currentIndex;
                neighbElementII.elementIndex[c5II] = neighbElementIII.currentIndex;

                // update the neighboures of neighbElementILeft and neighbElementIRight 
                neighbElementILeft.elementIndex[c7II] = neighbElementIRight.currentIndex;
                neighbElementIRight.elementIndex[c8II] = neighbElementILeft.currentIndex;                 

                // p2 now is eliminated.
                //posInList = quickSelect(pL, pI2, 0, pL.size()-1);
                //pL[posInList].keepPoint = false;
                pL[pI2].keepPoint = false;

                // eleminating the elements: currentElemet and neighbElementI
                currentElemet.keepElement = false;
                neighbElementI.keepElement = false;

                // updatong the all elOwner
                //posInList = quickSelect(pL, pI1, 0, pL.size()-1);
                //pL[posInList].elOwner = neighbElementILeft.currentIndex;
                //pL[pI3].elOwner = neighbElementIII.currentIndex;
                //pL[pI1].elOwner = neighbElementIII.currentIndex;
                //pL[pI4].elOwner = neighbElementIRight.currentIndex;
            }
        }//3
                 
    }//2

   Info << "                 After coarsening " << numberOfELDeleted << " element/elemens were deleted." << endl;
}//1

*/


/*
//- note: func to coarse the surface grid
template<class CloudType>
inline void Foam::FrontDynamic<CloudType>::frontCoarsening(DynamicList<pointData>& pL, DynamicList<elementInfo>& eL,
                                                           scalar minEdge, scalar maxEdge, scalar maxAspectRatio)
{//1
    label pI1, pI2, pI3, pI4;
    label c4I;
    //label c4II;
    label c4III;
    //label c5I
    label c5II;
    //label c5III;
    label c6I;
    label c6II;
    label c6III;
    //label c7I;
    label c7II;
    //label c7III;
    label c8I;
    label c8II;
    label c8III;
    label I, II, III;
    label f, cNI, cNII, cNIII;

    label numberOfELDeleted = 0;

    label initialElListSize = eL.size();

    //forAll(eL, elementI)
    for(int elementI=0; elementI < initialElListSize; elementI++)   
    //for(int elementI=0; elementI < initial_; elementI++)//JJ   
    {//2

        elementInfo& currentElemet = eL[elementI];

        // order the lengths, s1 shortest, s3 longest 
        scalar s1 = mag( currentElemet.points[1] - currentElemet.points[0]);
        scalar s2 = mag( currentElemet.points[2] - currentElemet.points[1]);
        scalar s3 = mag( currentElemet.points[0] - currentElemet.points[2]);
        label n1 = 0;
        label n2 = 1;
        label n3 = 2;
        bool smallElement = isElementSmall(currentElemet, minEdge, maxEdge, maxAspectRatio, s1, s2, s3, n1, n2, n3);
        Info<< "after smallElement " << ' ' << smallElement << ' ' << currentElemet.keepElement << endl;

        if ( smallElement && currentElemet.keepElement)
        {//3
            //Info<< "The element " << ' ' << currentElemet.currentIndex << ' ' << " will be coarsed." << endl;

            numberOfELDeleted = numberOfELDeleted + 2;

            label shortEdge = n1;

            // I, II, III are the point indexes of current element.
            // n1 is the index of shortest edge and p1.
            findPointsOrder( shortEdge, I, II, III);

            point p1 = currentElemet.points[I];
            point p2 = currentElemet.points[II];
            point p3 = currentElemet.points[III];
            pI1 = currentElemet.pointIndex[I];
            pI2 = currentElemet.pointIndex[II];
            pI3 = currentElemet.pointIndex[III];

            // find point 5, third point of neighbElementII
            //label posInList = quickSelect(eL, currentElemet.elementIndex[II], 0, eL.size()-1);
            //elementInfo& neighbElementII = eL[posInList];
            elementInfo& neighbElementII = eL[currentElemet.elementIndex[II]];
            for( f=0; f<=3; f++)
            {
                if( neighbElementII.pointIndex[f] != pI2 && neighbElementII.pointIndex[f] != pI3 )
                {
                    Info<< "f5 " << ' ' << f<< endl;
                    break;
                }
                if( f == 3)
                {
		    FatalErrorIn
		    (
			    " ^^^^^^^^^^^^^ we can not find such a point."
		    ) << abort(FatalError);
                }
            }
            // c5I, c5II, c5III are the point indexes of neighbElementII.
            findPointsOrder( f, cNI, cNII, cNIII);
            point p5 = neighbElementII.points[f];
            //c5I = cNI;
	    c5II = cNII;
	    //c5III = cNIII;

            // find point 6, third point of neighbElementIII
            //posInList = quickSelect(eL, currentElemet.elementIndex[III], 0, eL.size()-1);
            //elementInfo& neighbElementIII = eL[posInList];
            elementInfo& neighbElementIII = eL[currentElemet.elementIndex[III]];
            for( f=0; f<=3; f++)
            {
                if( neighbElementIII.pointIndex[f] != pI1 && neighbElementIII.pointIndex[f] != pI3 )
                {
                    Info<< "f6 " << ' ' << f<< endl;
                    break;
                }
                if( f == 3)
                {
		    FatalErrorIn
		    (
			    " ^^^^^^^^^^^^^ we can not find such a point."
		    ) << abort(FatalError);
                }
            }
            findPointsOrder( f, cNI, cNII, cNIII);
            point p6 = neighbElementIII.points[f];
            c6I = cNI;
	    c6II = cNII;
	    c6III = cNIII;

            // find point 4, third point of neighbElementIII
            //posInList = quickSelect(eL, currentElemet.elementIndex[I], 0, eL.size()-1);
            //elementInfo& neighbElementI = eL[posInList];
            elementInfo& neighbElementI = eL[currentElemet.elementIndex[I]];
            for( f=0; f<=3; f++)
            {
                if( neighbElementI.pointIndex[f] != pI1 && neighbElementI.pointIndex[f] != pI2 )
                {
                    Info<< "f4 " << ' ' << f<< endl;
                    break;
                }
                if( f == 3)
                {
		    FatalErrorIn
		    (
			    " ^^^^^^^^^^^^^ we can not find such a point."
		    ) << abort(FatalError);
                }
            }
            findPointsOrder( f, cNI, cNII, cNIII);
            point p4 = neighbElementI.points[f];
            pI4 = neighbElementI.pointIndex[f];
            c4I = cNI;
	    //c4II = cNII;
	    c4III = cNIII;

            // find point 7, third point of neighbElementILeft
            //posInList = quickSelect(eL, neighbElementI.elementIndex[c4III], 0, eL.size()-1);
            //elementInfo& neighbElementILeft = eL[posInList];
            elementInfo& neighbElementILeft = eL[neighbElementI.elementIndex[c4III]];
            for( f=0; f<=3; f++)
            {
                if( neighbElementILeft.pointIndex[f] != pI1 && neighbElementILeft.pointIndex[f] != pI4 )
                {
                    Info<< "f7 " << ' ' << f<< endl;
                    break;
                }
                if( f == 3)
                {
		    FatalErrorIn
		    (
			    " ^^^^^^^^^^^^^ we can not find such a point."
		    ) << abort(FatalError);
                }
            }
            findPointsOrder( f, cNI, cNII, cNIII);
            point p7 = neighbElementILeft.points[f];
            //c7I = cNI;
	    c7II = cNII;
	    //c7III = cNIII;

            // find point 8, third point of neighbElementIRight
            //posInList = quickSelect(eL, neighbElementI.elementIndex[c4I], 0, eL.size()-1);
            //elementInfo& neighbElementIRight = eL[posInList];
            elementInfo& neighbElementIRight = eL[neighbElementI.elementIndex[c4I]];
            for( f=0; f<=3; f++)
            {
                if( neighbElementIRight.pointIndex[f] != pI2 && neighbElementIRight.pointIndex[f] != pI4 )
                {
                    Info<< "f8 " << ' ' << f<< endl;
                    break;
                }
                if( f == 3)
                {
		    FatalErrorIn
		    (
			    " ^^^^^^^^^^^^^ we can not find such a point."
		    ) << abort(FatalError);
                }
            }
            findPointsOrder( f, cNI, cNII, cNIII);
            point p8 = neighbElementIRight.points[f];
            c8I = cNI;
	    c8II = cNII;
	    c8III = cNIII;

            // allocation of the new point: 
            point newPoint = 0.5 * ( p1 + p2) + 0.125 * ( p3 + p4) - 0.0625 * ( p5 + p6 + p7 + p8);

            // position of p1 now is updated.
            //posInList = quickSelect(pL, pI1, 0, pL.size()-1);
            //pL[posInList].currentPoint = newPoint;
             pL[pI1].currentPoint = newPoint;


            // transform only all positions of points of the element sharing the lable of p1 to newPoint 
            cNI = c6I;
            cNII = c6II;
            cNIII = c6III;              
            elementInfo elTemp1 = neighbElementIII;
            neighbElementIII.points[cNII] = newPoint;
            do 
            {
                //label posInList = quickSelect(eL, elTemp1.elementIndex[cNI], 0, eL.size()-1);
                //elementInfo& neighbElTemp1 = eL[posInList];
                elementInfo& neighbElTemp1 = eL[elTemp1.elementIndex[cNI]];
                for( f=0; f<=3; f++)
                {
                    if( neighbElTemp1.pointIndex[f] != elTemp1.pointIndex[cNI] && 
                        neighbElTemp1.pointIndex[f] != elTemp1.pointIndex[cNII] )
                    {
                        Info<< "F I can break " << ' ' << f << endl;
                        break;
                    }
                   if( f == 3)
                    {
		        FatalErrorIn
		        (
			    " ^^^^^^^^^^^^^ we can not find such a point."
		        ) << abort(FatalError);
                    }
                }
                findPointsOrder( f, cNI, cNII, cNIII);
                         
                elTemp1 = neighbElTemp1;
                neighbElTemp1.points[cNII] = newPoint;
                         //Info<< "F elTemp1.currentIndex/neighbElementILeft.currentIndex " << elTemp1.currentIndex << ' ' << neighbElementILeft.currentIndex << endl;
            } while( elTemp1.currentIndex != neighbElementILeft.currentIndex );
            //Info<< "First do " << ' ' << currentElemet.currentIndex << endl;

            // transform all positions and indexes of points of the element  sharing the p2 to newPoint
            // then the p2 must be eliminated.
            cNI = c8I;
            cNII = c8II;
            cNIII = c8III;              
            elementInfo elTemp2 = neighbElementIRight;
            neighbElementIRight.points[cNII] = newPoint;
            neighbElementIRight.pointIndex[cNII] = pI1;

                         Info<< "start/neighbour/currentElemet.currentIndex  " << elTemp2.currentIndex << ' '
                             << elTemp2.elementIndex[cNI] << ' ' << currentElemet.currentIndex << endl;
            Info<< "el/II/III/I/left/right  " << currentElemet.currentIndex << ' '
                << currentElemet.elementIndex[II] << ' ' << currentElemet.elementIndex[III] << ' ' << currentElemet.elementIndex[I] << ' ' 
                << neighbElementI.elementIndex[c4III] << ' ' << neighbElementI.elementIndex[c4I] << endl;
            do 
            {
                //label posInList = quickSelect(eL, elTemp2.elementIndex[cNI], 0, eL.size()-1);
                //elementInfo& neighbElTemp2 = eL[posInList];
                elementInfo& neighbElTemp2 = eL[elTemp2.elementIndex[cNI]];
                for( f=0; f<=3; f++)
                {
                    if( neighbElTemp2.pointIndex[f] != elTemp2.pointIndex[cNI] && 
                        neighbElTemp2.pointIndex[f] != elTemp2.pointIndex[cNII] )
                    {
                        Info<< "S I can break " << ' ' << f << endl;
                        break;
                    }
                    if( f == 3)
                    {
		        FatalErrorIn
		        (
			    " ^^^^^^^^^^^^^ we can not find such a point."
		        ) << abort(FatalError);
                    }
                }
                findPointsOrder( f, cNI, cNII, cNIII);

                elTemp2 = neighbElTemp2;
                neighbElTemp2.points[cNII] = newPoint;
                neighbElTemp2.pointIndex[cNII] = pI1;
                         Info<< "S elTemp2.currentIndex/neighbElementII.currentIndex " << elTemp2.currentIndex << ' ' << neighbElementII.currentIndex << endl;
            } while( elTemp2.currentIndex != neighbElementII.currentIndex );
            //Info<< "Second do " << ' ' << currentElemet.currentIndex << endl;
                     
            // p2 now is eliminated.
            //posInList = quickSelect(pL, pI2, 0, pL.size()-1);
            //pL[posInList].keepPoint = false;
            pL[pI2].keepPoint = false;
                  
            // update the neighboures of neighbElementIII and neighbElementII 
            neighbElementIII.elementIndex[c6II] = neighbElementII.currentIndex;
            neighbElementII.elementIndex[c5II] = neighbElementIII.currentIndex;

            // update the neighboures of neighbElementILeft and neighbElementIRight 
            neighbElementILeft.elementIndex[c7II] = neighbElementIRight.currentIndex;
            neighbElementIRight.elementIndex[c8II] = neighbElementILeft.currentIndex;                 

            // eleminating the elements: currentElemet and neighbElementI
            currentElemet.keepElement = false;
            neighbElementI.keepElement = false;

            // updatong the all elOwner
            //posInList = quickSelect(pL, pI1, 0, pL.size()-1);
            //pL[posInList].elOwner = neighbElementILeft.currentIndex;
            //pL[pI3].elOwner = neighbElementIII.currentIndex;
            //pL[pI1].elOwner = neighbElementIII.currentIndex;
            //pL[pI4].elOwner = neighbElementIRight.currentIndex;
        }//3
                 
    }//2

   Info << "                 After coarsening " << numberOfELDeleted << " element/elemens were deleted." << endl;
}//1
*/

//-note:  func to refine the surface grid
template<class CloudType>
inline void Foam::FrontDynamic<CloudType>::frontRefining( DynamicList<pointData>& pL, DynamicList<elementInfo>& eL,
                            scalar minEdge, scalar maxEdge, scalar maxAspectRatio)
{//1
    label pI1, pI2, pI3, pI4;
    label c4I;
    label c4II;
    label c4III;
    label c5I;
    label c5II;
    label c5III;
    label c6I;
    label c6II;
    label c6III;
    label c7I;
    label c7II;
    label c7III;
    label c8I;
    label c8II;
    label c8III;
    label c9I, c9II, c9III;
    label c10I, c10II, c10III;
    label c11I, c11II, c11III;
    label I, II, III;
    label f, cNI, cNII, cNIII;

    label numberOfELAdded = 0;
    
    label initialElListSize = eL.size();

    //forAll(eL, elementI)
    for(int elementI=0; elementI < initialElListSize; elementI++)      
    {//2
        elementInfo& currentElemet = eL[elementI];

            //Info<< "AHey " << ' ' << elementI << endl;

        // order the lengths, s1 shortest, s3 longest 
        scalar s1 = mag( currentElemet.points[1] - currentElemet.points[0]);
        scalar s2 = mag( currentElemet.points[2] - currentElemet.points[1]);
        scalar s3 = mag( currentElemet.points[0] - currentElemet.points[2]);
        label n1 = 0;
        label n2 = 1;
        label n3 = 2;
            //Info<< "s1/s2/s3 " << ' ' << s1 << ' ' << s2 << ' ' << s3 << endl;

        bool refinedElement = isElementRefined(currentElemet, minEdge, maxEdge, maxAspectRatio, s1, s2, s3, n1, n2, n3);

            //Info<< "AAHey--->large " << ' ' << refinedElement << ' ' << currentElemet.keepElement << ' ' << n3 << endl;

        if ( refinedElement && currentElemet.keepElement )
        {//3
            //Info<< "The element " << ' ' << currentElemet.currentIndex << ' ' << " will be refined." << endl;

            //pL.reserve(10*eL.size());
            //eL.reserve(10*eL.size());
             //Info<< "bef rserving " << endl;
            //pL.reserve( pL.size()+1);
            //eL.reserve( eL.size()+2);

            //numberOfELAdded = numberOfELAdded + 2;
             
            label largeEdge = n3;
                     
            // I, II, III are the point indexes of current element.
            // n1 is the index of shortest edge and p1.
            findPointsOrder( largeEdge, I, II, III);
            point p1 = currentElemet.points[I];//JJ
            point p2 = currentElemet.points[II];//JJ
            point p3 = currentElemet.points[III];//JJ
            pI1 = currentElemet.pointIndex[I];
            pI2 = currentElemet.pointIndex[II];
            pI3 = currentElemet.pointIndex[III];

            //Info << "all Nbs---> " << currentElemet.elementIndex[I] << ' ' << currentElemet.elementIndex[II] << ' ' 
            //     << currentElemet.elementIndex[III] << endl;

            // find point 5, third point of neighbElementII
            //label posInList = quickSelect(eL, currentElemet.elementIndex[II], 0, eL.size()-1);
            //elementInfo& neighbElementII = eL[posInList];
            elementInfo& neighbElementII = eL[currentElemet.elementIndex[II]];
            for( f=0; f<=3; f++)
            {
                if( neighbElementII.pointIndex[f] == pI2)
                {
                    //Info<< "f5 " << ' ' << f<< endl;
                    break;
                }
                if( f == 3)
                {
		    FatalErrorIn
		    (
			    " ^^^^^^^^^^^^^ we can not find such a point."
		    ) << abort(FatalError);
                }
            }            
            // c5I, c5II, c5III are the point indexes of neighbElementII.
            findPointsOrder( f, cNI, cNII, cNIII);
            point p5 = neighbElementII.points[cNII];
            c5I = cNII;
	    c5II = cNIII;
	    c5III = cNI;

            // find point 6, third point of neighbElementIII
            //posInList = quickSelect(eL, currentElemet.elementIndex[III], 0, eL.size()-1);
            //elementInfo& neighbElementIII = eL[posInList];
            elementInfo& neighbElementIII = eL[currentElemet.elementIndex[III]];
            for( f=0; f<=3; f++)
            {
                if( neighbElementIII.pointIndex[f] == pI3)
                {
                    //Info<< "f6 " << ' ' << f<< endl;
                    break;
                }
                if( f == 3)
                {
		    FatalErrorIn
		    (
			    " ^^^^^^^^^^^^^ we can not find such a point."
		    ) << abort(FatalError);
                }
            }
                    //Info<< "f6 " << ' ' << f<< endl;
            findPointsOrder( f, cNI, cNII, cNIII);
            point p6 = neighbElementIII.points[cNII];
	    c6I = cNII;
	    c6II = cNIII;
	    c6III = cNI;

            // find point 4, third point of neighbElementIII
            //posInList = quickSelect(eL, currentElemet.elementIndex[I], 0, eL.size()-1);
            //elementInfo& neighbElementI = eL[posInList];
            elementInfo& neighbElementI = eL[currentElemet.elementIndex[I]];
            for( f=0; f<=3; f++)
            {
                if( neighbElementI.pointIndex[f] == pI1)
                {
                    //Info<< "f4 " << ' ' << f<< endl;
                    break;
                }
                if( f == 3)
                {
		    FatalErrorIn
		    (
			    " ^^^^^^^^^^^^^ we can not find such a point."
		    ) << abort(FatalError);
                }
            }
            findPointsOrder( f, cNI, cNII, cNIII);
            point p4 = neighbElementI.points[cNII];
            pI4 = neighbElementI.pointIndex[cNII];
	    c4I = cNII;
	    c4II = cNIII;
	    c4III = cNI;


            // find point 7, third point of neighbElementILeft
            //posInList = quickSelect(eL, neighbElementI.elementIndex[c4III], 0, eL.size()-1);
            //elementInfo& neighbElementILeft = eL[posInList];
            elementInfo& neighbElementILeft = eL[neighbElementI.elementIndex[c4III]];
            for( f=0; f<=3; f++)
            {
                if( neighbElementILeft.pointIndex[f] == pI1)
                {
                    //Info<< "f7 " << ' ' << f<< endl;
                    break;
                }
                if( f == 3)
                {
		    FatalErrorIn
		    (
			    " ^^^^^^^^^^^^^ we can not find such a point."
		    ) << abort(FatalError);
                }
            }
            findPointsOrder( f, cNI, cNII, cNIII);
            point p7 = neighbElementILeft.points[cNII];
	    c7I = cNII;
	    c7II = cNIII;
	    c7III = cNI;
                    //Info<< "after f7 " << ' ' << f<< endl;


                         //Info << "el/II/III/I/left/right  " << currentElemet.currentIndex << ' '
                         //     << currentElemet.elementIndex[II] << ' ' << currentElemet.elementIndex[III] << ' ' << currentElemet.elementIndex[I] << ' ' 
                          //    << neighbElementI.elementIndex[c4III] << ' ' << neighbElementI.elementIndex[c4I] << endl;

            // find point 8, third point of neighbElementIRight
            //posInList = quickSelect(eL, neighbElementI.elementIndex[c4I], 0, eL.size()-1);
            //elementInfo& neighbElementIRight = eL[posInList];
            elementInfo& neighbElementIRight = eL[neighbElementI.elementIndex[c4I]];
                    //Info<< "neighbElementIRight " << ' ' << neighbElementI.elementIndex[c4I] << endl;

            for( f=0; f<=3; f++)
            {
                if( neighbElementIRight.pointIndex[f] == pI4)
                {
                    //Info<< "f8 " << ' ' << f<< endl;
                    break;
                }
                if( f == 3)
                {
		    FatalErrorIn
		    (
			    " ^^^^^^^^^^^^^ we can not find such a point."
		    ) << abort(FatalError);
                }
            }
            findPointsOrder( f, cNI, cNII, cNIII);
            point p8 = neighbElementIRight.points[cNII];
	    c8I = cNII;
	    c8II = cNIII;
	    c8III = cNI;

                    numberOfELAdded = numberOfELAdded + 2;

		    // allocation of the new point:
		    pointData newPoint;
		    newPoint.currentIndex = pL.size();//pL[pL.size()-1].currentIndex + 1
		    newPoint.currentPoint = 0.5 * ( p1 + p2) + 0.125 * ( p3 + p4) - 0.0625 * ( p5 + p6 + p7 + p8);
		    newPoint.keepPoint = true;
		    newPoint.elOwner = currentElemet.currentIndex;

		    //Info<< "   p1/p2/p3/p4/p5/p6/p7/p8: " << p1 << ' ' << p2 << ' ' << p3 << ' ' << p4 << ' '
		    //                                      << p5 << ' ' << p6 << ' ' << p7 << ' ' << p8 << endl;
		    //Info<< "newPoint for refining " << ' ' << newPoint.currentPoint << endl;

		    // allocation of the new element at the upper part:
		    elementInfo newElementU;
		    newElementU.currentIndex = eL.size();//eL[eL.size()-1].currentIndex + 1;
		    newElementU.points[0] = newPoint.currentPoint;
		    newElementU.points[1] = p3;
		    newElementU.points[2] = p1;
		    newElementU.pointIndex[0] = newPoint.currentIndex;
		    newElementU.pointIndex[1] = pI3;
		    newElementU.pointIndex[2] = pI1;
		    newElementU.keepElement = true;

		    //Info<< "after newElementU " << ' ' << eL.size() << endl;

		    // allocation of the new element at the bottom part:
		    elementInfo newElementB;
		    newElementB.currentIndex = eL.size()+1;//eL[eL.size()-1].currentIndex + 1;
		    newElementB.points[0] = newPoint.currentPoint;
		    newElementB.points[1] = p1;
		    newElementB.points[2] = p4;
		    newElementB.pointIndex[0] = newPoint.currentIndex;
		    newElementB.pointIndex[1] = pI1;
		    newElementB.pointIndex[2] = pI4;
		    newElementB.keepElement = true;


		    //Info<< "after newElementB " << ' ' << eL.size() << endl;

		    // transform currentElemet[I] and neighbElementI[c4III] to point to the newPoint 
		    currentElemet.points[I] = newPoint.currentPoint;
		    currentElemet.pointIndex[I] = newPoint.currentIndex;
		    neighbElementI.points[c4III] = newPoint.currentPoint;
		    neighbElementI.pointIndex[c4III] = newPoint.currentIndex;

		    //Info<< "after transform currentElemet[I]  " << ' ' << eL.size() << endl;

		    // update the neighboures of newElementU
		    newElementU.elementIndex[0] = currentElemet.currentIndex;
		    newElementU.elementIndex[1] = neighbElementIII.currentIndex;
		    newElementU.elementIndex[2] = newElementB.currentIndex;

		    //Info<< "after update the neighboures of newElementU " << ' ' << eL.size() << endl;

		    // update the neighboures of newElementB
		    newElementB.elementIndex[0] = newElementU.currentIndex;
		    newElementB.elementIndex[1] = neighbElementILeft.currentIndex;
		    newElementB.elementIndex[2] = neighbElementI.currentIndex;

		    //Info<< "after update the neighboures of newElementB " << ' ' << eL.size() << endl;

		    // update the neighboures of currentElemet and neighbElementI
		    currentElemet.elementIndex[III] = newElementU.currentIndex;
		    neighbElementI.elementIndex[c4III] = newElementB.currentIndex;

		    //Info<< "after update the neighboures of currentElemet and neighbElementI " << ' ' << eL.size() << endl;

		    // update the neighboures of neighbElementIII and neighbElementILeft
		    neighbElementIII.elementIndex[c6II] = newElementU.currentIndex;
		    neighbElementILeft.elementIndex[c7II] = newElementB.currentIndex;

		    //Info<< "after update the neighboures of neighbElementIII and neighbElementILeft " << ' ' << eL.size() << endl;

		    //updating the elOwner for p1
		    //posInList = quickSelect(pL, pI1, 0, pL.size()-1);
		    //pL[posInList].elOwner = neighbElementIII.currentIndex;
		    pL[pI1].elOwner = neighbElementILeft.currentIndex;

		    //Info << "--$$$---@@@@@28@@@@----neighbElementIII.elementIndex[c6I] " << neighbElementIII.elementIndex[0] << ' '
		    //     << neighbElementIII.currentIndex << endl;
		    //Info << "--$$$---@@@@@28@@@@----neighbElementILeft.elementIndex[c7I] " << neighbElementILeft.elementIndex[0] << ' '
		    //     << neighbElementILeft.currentIndex << endl;

		   // Info << "--$$$---@@@@@28@@@@----neighbElementIII.elementIndex[c6II] " << neighbElementIII.elementIndex[1] << ' '
		   //      << neighbElementIII.currentIndex << endl;
		   // Info << "--$$$---@@@@@28@@@@----neighbElementILeft.elementIndex[c7II] " << neighbElementILeft.elementIndex[1] << ' '
		   //      << neighbElementILeft.currentIndex << endl;

		   // Info << "--$$$---@@@@@28@@@@----neighbElementIII.elementIndex[c6III] " << neighbElementIII.elementIndex[2] << ' '
		   //      << neighbElementIII.currentIndex << endl;
		   // Info << "--$$$---@@@@@28@@@@----neighbElementILeft.elementIndex[c7III] " << neighbElementILeft.elementIndex[2] << ' '
		   //      << neighbElementILeft.currentIndex << endl;

		    pL.append(newPoint);
		    eL.append(newElementU);
		    eL.append(newElementB);
            //Info<< "FFFFFFFFFFFffinal " << ' ' << elementI << endl;
       }//3
   }//2
   Info << "                 After refining " << numberOfELAdded << " element/elemens were added." << endl;
}//1

/*
//-note:  func to refine the surface grid
template<class CloudType>
inline void Foam::FrontDynamic<CloudType>::frontRefining( DynamicList<pointData>& pL, DynamicList<elementInfo>& eL,
                            scalar minEdge, scalar maxEdge, scalar maxAspectRatio)
{//1
    label pI1, pI2, pI3, pI4;
    label c4I;
    label c4II;
    label c4III;
    label c5I;
    label c5II;
    label c5III;
    label c6I;
    label c6II;
    label c6III;
    label c7I;
    label c7II;
    label c7III;
    label c8I;
    label c8II;
    label c8III;
    label c9I, c9II, c9III;
    label c10I, c10II, c10III;
    label c11I, c11II, c11III;
    label I, II, III;
    label f, cNI, cNII, cNIII;

    label numberOfELAdded = 0;
    
    label initialElListSize = eL.size();

    //forAll(eL, elementI)
    for(int elementI=0; elementI < initialElListSize; elementI++)      
    {//2
        elementInfo& currentElemet = eL[elementI];

            //Info<< "AHey " << ' ' << elementI << endl;

        // order the lengths, s1 shortest, s3 longest 
        scalar s1 = mag( currentElemet.points[1] - currentElemet.points[0]);
        scalar s2 = mag( currentElemet.points[2] - currentElemet.points[1]);
        scalar s3 = mag( currentElemet.points[0] - currentElemet.points[2]);
        label n1 = 0;
        label n2 = 1;
        label n3 = 2;
            //Info<< "s1/s2/s3 " << ' ' << s1 << ' ' << s2 << ' ' << s3 << endl;

        bool refinedElement = isElementRefined(currentElemet, minEdge, maxEdge, maxAspectRatio, s1, s2, s3, n1, n2, n3);

            //Info<< "AAHey--->large " << ' ' << refinedElement << ' ' << currentElemet.keepElement << ' ' << n3 << endl;

        if ( refinedElement && currentElemet.keepElement )
        {//3
            //Info<< "The element " << ' ' << currentElemet.currentIndex << ' ' << " will be refined." << endl;

            //pL.reserve(10*eL.size());
            //eL.reserve(10*eL.size());
             //Info<< "bef rserving " << endl;
            //pL.reserve( pL.size()+1);
            //eL.reserve( eL.size()+2);

            //numberOfELAdded = numberOfELAdded + 2;
             
            label largeEdge = n3;
                     
            // I, II, III are the point indexes of current element.
            // n1 is the index of shortest edge and p1.
            findPointsOrder( largeEdge, I, II, III);
            point p1 = currentElemet.points[I];//JJ
            point p2 = currentElemet.points[II];//JJ
            point p3 = currentElemet.points[III];//JJ
            pI1 = currentElemet.pointIndex[I];
            pI2 = currentElemet.pointIndex[II];
            pI3 = currentElemet.pointIndex[III];

            //Info << "all Nbs---> " << currentElemet.elementIndex[I] << ' ' << currentElemet.elementIndex[II] << ' ' 
            //     << currentElemet.elementIndex[III] << endl;

            // find point 5, third point of neighbElementII
            //label posInList = quickSelect(eL, currentElemet.elementIndex[II], 0, eL.size()-1);
            //elementInfo& neighbElementII = eL[posInList];
            elementInfo& neighbElementII = eL[currentElemet.elementIndex[II]];
            for( f=0; f<=3; f++)
            {
                if( neighbElementII.pointIndex[f] == pI2)
                {
                    //Info<< "f5 " << ' ' << f<< endl;
                    break;
                }
                if( f == 3)
                {
		    FatalErrorIn
		    (
			    " ^^^^^^^^^^^^^ we can not find such a point."
		    ) << abort(FatalError);
                }
            }            
            // c5I, c5II, c5III are the point indexes of neighbElementII.
            findPointsOrder( f, cNI, cNII, cNIII);
            point p5 = neighbElementII.points[cNII];
            c5I = cNII;
	    c5II = cNIII;
	    c5III = cNI;

            // find point 6, third point of neighbElementIII
            //posInList = quickSelect(eL, currentElemet.elementIndex[III], 0, eL.size()-1);
            //elementInfo& neighbElementIII = eL[posInList];
            elementInfo& neighbElementIII = eL[currentElemet.elementIndex[III]];
            for( f=0; f<=3; f++)
            {
                if( neighbElementIII.pointIndex[f] == pI3)
                {
                    //Info<< "f6 " << ' ' << f<< endl;
                    break;
                }
                if( f == 3)
                {
		    FatalErrorIn
		    (
			    " ^^^^^^^^^^^^^ we can not find such a point."
		    ) << abort(FatalError);
                }
            }
                    //Info<< "f6 " << ' ' << f<< endl;
            findPointsOrder( f, cNI, cNII, cNIII);
            point p6 = neighbElementIII.points[cNII];
	    c6I = cNII;
	    c6II = cNIII;
	    c6III = cNI;

            // find point 4, third point of neighbElementIII
            //posInList = quickSelect(eL, currentElemet.elementIndex[I], 0, eL.size()-1);
            //elementInfo& neighbElementI = eL[posInList];
            elementInfo& neighbElementI = eL[currentElemet.elementIndex[I]];
            for( f=0; f<=3; f++)
            {
                if( neighbElementI.pointIndex[f] == pI1)
                {
                    //Info<< "f4 " << ' ' << f<< endl;
                    break;
                }
                if( f == 3)
                {
		    FatalErrorIn
		    (
			    " ^^^^^^^^^^^^^ we can not find such a point."
		    ) << abort(FatalError);
                }
            }
            findPointsOrder( f, cNI, cNII, cNIII);
            point p4 = neighbElementI.points[cNII];
            pI4 = neighbElementI.pointIndex[cNII];
	    c4I = cNII;
	    c4II = cNIII;
	    c4III = cNI;


            // find point 7, third point of neighbElementILeft
            //posInList = quickSelect(eL, neighbElementI.elementIndex[c4III], 0, eL.size()-1);
            //elementInfo& neighbElementILeft = eL[posInList];
            elementInfo& neighbElementILeft = eL[neighbElementI.elementIndex[c4III]];
            for( f=0; f<=3; f++)
            {
                if( neighbElementILeft.pointIndex[f] == pI1)
                {
                    //Info<< "f7 " << ' ' << f<< endl;
                    break;
                }
                if( f == 3)
                {
		    FatalErrorIn
		    (
			    " ^^^^^^^^^^^^^ we can not find such a point."
		    ) << abort(FatalError);
                }
            }
            findPointsOrder( f, cNI, cNII, cNIII);
            point p7 = neighbElementILeft.points[cNII];
	    c7I = cNII;
	    c7II = cNIII;
	    c7III = cNI;
                    //Info<< "after f7 " << ' ' << f<< endl;


                         //Info << "el/II/III/I/left/right  " << currentElemet.currentIndex << ' '
                         //     << currentElemet.elementIndex[II] << ' ' << currentElemet.elementIndex[III] << ' ' << currentElemet.elementIndex[I] << ' ' 
                          //    << neighbElementI.elementIndex[c4III] << ' ' << neighbElementI.elementIndex[c4I] << endl;

            // find point 8, third point of neighbElementIRight
            //posInList = quickSelect(eL, neighbElementI.elementIndex[c4I], 0, eL.size()-1);
            //elementInfo& neighbElementIRight = eL[posInList];
            elementInfo& neighbElementIRight = eL[neighbElementI.elementIndex[c4I]];
                    //Info<< "neighbElementIRight " << ' ' << neighbElementI.elementIndex[c4I] << endl;

            for( f=0; f<=3; f++)
            {
                if( neighbElementIRight.pointIndex[f] == pI4)
                {
                    //Info<< "f8 " << ' ' << f<< endl;
                    break;
                }
                if( f == 3)
                {
		    FatalErrorIn
		    (
			    " ^^^^^^^^^^^^^ we can not find such a point."
		    ) << abort(FatalError);
                }
            }
            findPointsOrder( f, cNI, cNII, cNIII);
            point p8 = neighbElementIRight.points[cNII];
	    c8I = cNII;
	    c8II = cNIII;
	    c8III = cNI;

            if( neighbElementIRight.currentIndex == neighbElementII.currentIndex  )
            {
                Info << " RRRRRRRRRRRRRRRRRRRRR---the right triple is occured." << currentElemet.currentIndex << endl;
                currentElemet.points[II] = p4;
                currentElemet.pointIndex[II] = pI4;

                // finding neighbElementIIRight
                elementInfo& neighbElementIIRight = eL[neighbElementII.elementIndex[c5I]];
                for( f=0; f<=3; f++)
                {
                    if( neighbElementIIRight.pointIndex[f] == pI4)
                    {
                        break;
                    }
                    if( f == 3)
                    {
		        FatalErrorIn
		        (
			    " ^^^^^^^^^^^^^ we can not find pI4 for neighbElementIIRight."
		        ) << abort(FatalError);
                    }
                }
                findPointsOrder( f, cNI, cNII, cNIII);
	        c9I = cNII;
	        c9II = cNIII;
	        c9III = cNI;
                neighbElementIIRight.elementIndex[c9II] = currentElemet.currentIndex;
                currentElemet.elementIndex[II] = neighbElementIIRight.currentIndex;

                // finding neighbElementIBottom
                elementInfo& neighbElementIBottom = eL[neighbElementI.elementIndex[c4III]];
                for( f=0; f<=3; f++)
                {
                    if( neighbElementIBottom.pointIndex[f] == pI1)
                    {
                        break;
                    }
                    if( f == 3)
                    {
		        FatalErrorIn
		        (
			    " ^^^^^^^^^^^^^ we can not find pI1 for neighbElementIBottom."
		        ) << abort(FatalError);
                    }
                }
                findPointsOrder( f, cNI, cNII, cNIII);
	        c10I = cNII;
	        c10II = cNIII;
	        c10III = cNI;
                neighbElementIBottom.elementIndex[c10II] = currentElemet.currentIndex;
                currentElemet.elementIndex[I] = neighbElementIBottom.currentIndex;

                pL[pI2].keepPoint = false;
                neighbElementI.keepElement = false;
                neighbElementII.keepElement = false;
            }
            else if( neighbElementILeft.currentIndex == neighbElementIII.currentIndex  )
            {
                Info << " RRRRRRRRRRRRRRRRRRRRR---the left triple is occured." << currentElemet.currentIndex << endl;
                currentElemet.points[I] = p4;
                currentElemet.pointIndex[I] = pI4;

                // finding neighbElementIIILeft
                elementInfo& neighbElementIIILeft = eL[neighbElementIII.elementIndex[c6III]];
                for( f=0; f<=3; f++)
                {
                    if( neighbElementIIILeft.pointIndex[f] == pI3)
                    {
                        break;
                    }
                    if( f == 3)
                    {
		        FatalErrorIn
		        (
			    " ^^^^^^^^^^^^^ we can not find pI3 neighbElementIIILeft."
		        ) << abort(FatalError);
                    }
                }
                findPointsOrder( f, cNI, cNII, cNIII);
	        c11I = cNII;
	        c11II = cNIII;
	        c11III = cNI;
                neighbElementIIILeft.elementIndex[c11II] = currentElemet.currentIndex;
                currentElemet.elementIndex[III] = neighbElementIIILeft.currentIndex;

                // finding neighbElementIBottom
                elementInfo& neighbElementIBottom = eL[neighbElementI.elementIndex[c4I]];
                for( f=0; f<=3; f++)
                {
                    if( neighbElementIBottom.pointIndex[f] == pI4)
                    {
                        break;
                    }
                    if( f == 3)
                    {
		        FatalErrorIn
		        (
			    " ^^^^^^^^^^^^^ we can not find pI4 for neighbElementIBottom."
		        ) << abort(FatalError);
                    }
                }
                findPointsOrder( f, cNI, cNII, cNIII);
	        c10I = cNII;
	        c10II = cNIII;
	        c10III = cNI;
                neighbElementIBottom.elementIndex[c10II] = currentElemet.currentIndex;
                currentElemet.elementIndex[I] = neighbElementIBottom.currentIndex;

                pL[pI1].keepPoint = false;
                neighbElementI.keepElement = false;
                neighbElementIII.keepElement = false;
            }
            else
            {
                    numberOfELAdded = numberOfELAdded + 2;

		    // allocation of the new point:
		    pointData newPoint;
		    newPoint.currentIndex = pL.size();//pL[pL.size()-1].currentIndex + 1
		    newPoint.currentPoint = 0.5 * ( p1 + p2) + 0.125 * ( p3 + p4) - 0.0625 * ( p5 + p6 + p7 + p8);
		    newPoint.keepPoint = true;
		    newPoint.elOwner = currentElemet.currentIndex;

		    //Info<< "   p1/p2/p3/p4/p5/p6/p7/p8: " << p1 << ' ' << p2 << ' ' << p3 << ' ' << p4 << ' '
		    //                                      << p5 << ' ' << p6 << ' ' << p7 << ' ' << p8 << endl;
		    //Info<< "newPoint for refining " << ' ' << newPoint.currentPoint << endl;

		    // allocation of the new element at the upper part:
		    elementInfo newElementU;
		    newElementU.currentIndex = eL.size();//eL[eL.size()-1].currentIndex + 1;
		    newElementU.points[0] = newPoint.currentPoint;
		    newElementU.points[1] = p3;
		    newElementU.points[2] = p1;
		    newElementU.pointIndex[0] = newPoint.currentIndex;
		    newElementU.pointIndex[1] = pI3;
		    newElementU.pointIndex[2] = pI1;
		    newElementU.keepElement = true;

		    //Info<< "after newElementU " << ' ' << eL.size() << endl;

		    // allocation of the new element at the bottom part:
		    elementInfo newElementB;
		    newElementB.currentIndex = eL.size()+1;//eL[eL.size()-1].currentIndex + 1;
		    newElementB.points[0] = newPoint.currentPoint;
		    newElementB.points[1] = p1;
		    newElementB.points[2] = p4;
		    newElementB.pointIndex[0] = newPoint.currentIndex;
		    newElementB.pointIndex[1] = pI1;
		    newElementB.pointIndex[2] = pI4;
		    newElementB.keepElement = true;


		    //Info<< "after newElementB " << ' ' << eL.size() << endl;

		    // transform currentElemet[I] and neighbElementI[c4III] to point to the newPoint 
		    currentElemet.points[I] = newPoint.currentPoint;
		    currentElemet.pointIndex[I] = newPoint.currentIndex;
		    neighbElementI.points[c4III] = newPoint.currentPoint;
		    neighbElementI.pointIndex[c4III] = newPoint.currentIndex;

		    //Info<< "after transform currentElemet[I]  " << ' ' << eL.size() << endl;

		    // update the neighboures of newElementU
		    newElementU.elementIndex[0] = currentElemet.currentIndex;
		    newElementU.elementIndex[1] = neighbElementIII.currentIndex;
		    newElementU.elementIndex[2] = newElementB.currentIndex;

		    //Info<< "after update the neighboures of newElementU " << ' ' << eL.size() << endl;

		    // update the neighboures of newElementB
		    newElementB.elementIndex[0] = newElementU.currentIndex;
		    newElementB.elementIndex[1] = neighbElementILeft.currentIndex;
		    newElementB.elementIndex[2] = neighbElementI.currentIndex;

		    //Info<< "after update the neighboures of newElementB " << ' ' << eL.size() << endl;

		    // update the neighboures of currentElemet and neighbElementI
		    currentElemet.elementIndex[III] = newElementU.currentIndex;
		    neighbElementI.elementIndex[c4III] = newElementB.currentIndex;

		    //Info<< "after update the neighboures of currentElemet and neighbElementI " << ' ' << eL.size() << endl;

		    // update the neighboures of neighbElementIII and neighbElementILeft
		    neighbElementIII.elementIndex[c6II] = newElementU.currentIndex;
		    neighbElementILeft.elementIndex[c7II] = newElementB.currentIndex;

		    //Info<< "after update the neighboures of neighbElementIII and neighbElementILeft " << ' ' << eL.size() << endl;

		    //updating the elOwner for p1
		    //posInList = quickSelect(pL, pI1, 0, pL.size()-1);
		    //pL[posInList].elOwner = neighbElementIII.currentIndex;
		    pL[pI1].elOwner = neighbElementILeft.currentIndex;

		    //Info << "--$$$---@@@@@28@@@@----neighbElementIII.elementIndex[c6I] " << neighbElementIII.elementIndex[0] << ' '
		    //     << neighbElementIII.currentIndex << endl;
		    //Info << "--$$$---@@@@@28@@@@----neighbElementILeft.elementIndex[c7I] " << neighbElementILeft.elementIndex[0] << ' '
		    //     << neighbElementILeft.currentIndex << endl;

		   // Info << "--$$$---@@@@@28@@@@----neighbElementIII.elementIndex[c6II] " << neighbElementIII.elementIndex[1] << ' '
		   //      << neighbElementIII.currentIndex << endl;
		   // Info << "--$$$---@@@@@28@@@@----neighbElementILeft.elementIndex[c7II] " << neighbElementILeft.elementIndex[1] << ' '
		   //      << neighbElementILeft.currentIndex << endl;

		   // Info << "--$$$---@@@@@28@@@@----neighbElementIII.elementIndex[c6III] " << neighbElementIII.elementIndex[2] << ' '
		   //      << neighbElementIII.currentIndex << endl;
		   // Info << "--$$$---@@@@@28@@@@----neighbElementILeft.elementIndex[c7III] " << neighbElementILeft.elementIndex[2] << ' '
		   //      << neighbElementILeft.currentIndex << endl;

		    pL.append(newPoint);
		    eL.append(newElementU);
		    eL.append(newElementB);
            }
            //Info<< "FFFFFFFFFFFffinal " << ' ' << elementI << endl;
       }//3
   }//2

   Info << "                 After refining " << numberOfELAdded << " element/elemens were added." << endl;
}//1
*/

/*
//-note:  func to refine the surface grid
template<class CloudType>
inline void Foam::FrontDynamic<CloudType>::frontRefining( DynamicList<pointData>& pL, DynamicList<elementInfo>& eL,
                            scalar minEdge, scalar maxEdge, scalar maxAspectRatio)
{//1
    label pI1, pI2, pI3, pI4;
    label c4I;
    //label c4II;
    label c4III;
    //label c5I
    //label c5II;
    //label c5III;
    label c6I;
    label c6II;
    label c6III;
    label c7I;
    label c7II;
    label c7III;
    //label c8I;
    //label c8II;
    //label c8III;
    label I, II, III;
    label f, cNI, cNII, cNIII;

    label numberOfELAdded = 0;
    
    label initialElListSize = eL.size();

    //forAll(eL, elementI)
    for(int elementI=0; elementI < initialElListSize; elementI++)      
    {//2
        elementInfo& currentElemet = eL[elementI];

            Info<< "AHey " << ' ' << elementI << endl;

        // order the lengths, s1 shortest, s3 longest 
        scalar s1 = mag( currentElemet.points[1] - currentElemet.points[0]);
        scalar s2 = mag( currentElemet.points[2] - currentElemet.points[1]);
        scalar s3 = mag( currentElemet.points[0] - currentElemet.points[2]);
        label n1 = 0;
        label n2 = 1;
        label n3 = 2;
            Info<< "s1/s2/s3 " << ' ' << s1 << ' ' << s2 << ' ' << s3 << endl;

        bool refinedElement = isElementRefined(currentElemet, minEdge, maxEdge, maxAspectRatio, s1, s2, s3, n1, n2, n3);

            Info<< "AAHey--->large " << ' ' << refinedElement << ' ' << currentElemet.keepElement << ' ' << n3 << endl;

        if ( refinedElement && currentElemet.keepElement )
        {//3
            //Info<< "The element " << ' ' << currentElemet.currentIndex << ' ' << " will be refined." << endl;

            //pL.reserve(10*eL.size());
            //eL.reserve(10*eL.size());
             Info<< "bef rserving " << endl;
            //pL.reserve( pL.size()+1);
            //eL.reserve( eL.size()+2);

            numberOfELAdded = numberOfELAdded + 2;
             
            label largeEdge = n3;
                     
            // I, II, III are the point indexes of current element.
            // n1 is the index of shortest edge and p1.
            findPointsOrder( largeEdge, I, II, III);
            point p1 = currentElemet.points[I];//JJ
            point p2 = currentElemet.points[II];//JJ
            point p3 = currentElemet.points[III];//JJ
            pI1 = currentElemet.pointIndex[I];
            pI2 = currentElemet.pointIndex[II];
            pI3 = currentElemet.pointIndex[III];

            Info << "all Nbs---> " << currentElemet.elementIndex[I] << ' ' << currentElemet.elementIndex[II] << ' ' 
                 << currentElemet.elementIndex[III] << endl;

            // find point 5, third point of neighbElementII
            //label posInList = quickSelect(eL, currentElemet.elementIndex[II], 0, eL.size()-1);
            //elementInfo& neighbElementII = eL[posInList];
            elementInfo& neighbElementII = eL[currentElemet.elementIndex[II]];
            for( f=0; f<=3; f++)
            {
                if( neighbElementII.pointIndex[f] == pI2)
                {
                    Info<< "f5 " << ' ' << f<< endl;
                    break;
                }
                if( f == 3)
                {
		    FatalErrorIn
		    (
			    " ^^^^^^^^^^^^^ we can not find such a point."
		    ) << abort(FatalError);
                }
            }            
            // c5I, c5II, c5III are the point indexes of neighbElementII.
            findPointsOrder( f, cNI, cNII, cNIII);
            point p5 = neighbElementII.points[cNII];
            //c5I = cNI;
	    //c5II = cNII;
	    //c5III = cNIII;

            // find point 6, third point of neighbElementIII
            //posInList = quickSelect(eL, currentElemet.elementIndex[III], 0, eL.size()-1);
            //elementInfo& neighbElementIII = eL[posInList];
            elementInfo& neighbElementIII = eL[currentElemet.elementIndex[III]];
            for( f=0; f<=3; f++)
            {
                if( neighbElementIII.pointIndex[f] == pI3)
                {
                    Info<< "f6 " << ' ' << f<< endl;
                    break;
                }
                if( f == 3)
                {
		    FatalErrorIn
		    (
			    " ^^^^^^^^^^^^^ we can not find such a point."
		    ) << abort(FatalError);
                }
            }
                    Info<< "f6 " << ' ' << f<< endl;
            findPointsOrder( f, cNI, cNII, cNIII);
            point p6 = neighbElementIII.points[cNII];
	    c6I = cNII;
	    c6II = cNIII;
	    c6III = cNI;

            // find point 4, third point of neighbElementIII
            //posInList = quickSelect(eL, currentElemet.elementIndex[I], 0, eL.size()-1);
            //elementInfo& neighbElementI = eL[posInList];
            elementInfo& neighbElementI = eL[currentElemet.elementIndex[I]];
            for( f=0; f<=3; f++)
            {
                if( neighbElementI.pointIndex[f] == pI1)
                {
                    Info<< "f4 " << ' ' << f<< endl;
                    break;
                }
                if( f == 3)
                {
		    FatalErrorIn
		    (
			    " ^^^^^^^^^^^^^ we can not find such a point."
		    ) << abort(FatalError);
                }
            }
            findPointsOrder( f, cNI, cNII, cNIII);
            point p4 = neighbElementI.points[cNII];
            pI4 = neighbElementI.pointIndex[cNII];
	    c4I = cNII;
	    //c4II = cNII;
	    c4III = cNI;


            // find point 7, third point of neighbElementILeft
            //posInList = quickSelect(eL, neighbElementI.elementIndex[c4III], 0, eL.size()-1);
            //elementInfo& neighbElementILeft = eL[posInList];
            elementInfo& neighbElementILeft = eL[neighbElementI.elementIndex[c4III]];
            for( f=0; f<=3; f++)
            {
                if( neighbElementILeft.pointIndex[f] == pI1)
                {
                    Info<< "f7 " << ' ' << f<< endl;
                    break;
                }
                if( f == 3)
                {
		    FatalErrorIn
		    (
			    " ^^^^^^^^^^^^^ we can not find such a point."
		    ) << abort(FatalError);
                }
            }
            findPointsOrder( f, cNI, cNII, cNIII);
            point p7 = neighbElementILeft.points[cNII];
	    c7I = cNII;
	    c7II = cNIII;
	    c7III = cNI;
                    Info<< "after f7 " << ' ' << f<< endl;


                         Info << "el/II/III/I/left/right  " << currentElemet.currentIndex << ' '
                              << currentElemet.elementIndex[II] << ' ' << currentElemet.elementIndex[III] << ' ' << currentElemet.elementIndex[I] << ' ' 
                              << neighbElementI.elementIndex[c4III] << ' ' << neighbElementI.elementIndex[c4I] << endl;

            // find point 8, third point of neighbElementIRight
            //posInList = quickSelect(eL, neighbElementI.elementIndex[c4I], 0, eL.size()-1);
            //elementInfo& neighbElementIRight = eL[posInList];
            elementInfo& neighbElementIRight = eL[neighbElementI.elementIndex[c4I]];
                    Info<< "neighbElementIRight " << ' ' << neighbElementI.elementIndex[c4I] << endl;

            for( f=0; f<=3; f++)
            {
                if( neighbElementIRight.pointIndex[f] == pI4)
                {
                    Info<< "f8 " << ' ' << f<< endl;
                    break;
                }
                if( f == 3)
                {
		    FatalErrorIn
		    (
			    " ^^^^^^^^^^^^^ we can not find such a point."
		    ) << abort(FatalError);
                }
            }
            findPointsOrder( f, cNI, cNII, cNIII);
            point p8 = neighbElementIRight.points[cNII];
	    //c8I = cNI;
	    //c8II = cNII;
	    //c8III = cNIII;

            // allocation of the new point:
            pointData newPoint;
            newPoint.currentIndex = pL.size();//pL[pL.size()-1].currentIndex + 1
            newPoint.currentPoint = 0.5 * ( p1 + p2) + 0.125 * ( p3 + p4) - 0.0625 * ( p5 + p6 + p7 + p8);
            newPoint.keepPoint = true;
            newPoint.elOwner = currentElemet.currentIndex;

            Info<< "   p1/p2/p3/p4/p5/p6/p7/p8: " << p1 << ' ' << p2 << ' ' << p3 << ' ' << p4 << ' '
                                                  << p5 << ' ' << p6 << ' ' << p7 << ' ' << p8 << endl;
            Info<< "newPoint for refining " << ' ' << newPoint.currentPoint << endl;

            // allocation of the new element at the upper part:
            elementInfo newElementU;
            newElementU.currentIndex = eL.size();//eL[eL.size()-1].currentIndex + 1;
            newElementU.points[0] = newPoint.currentPoint;
            newElementU.points[1] = p3;
            newElementU.points[2] = p1;
            newElementU.pointIndex[0] = newPoint.currentIndex;
            newElementU.pointIndex[1] = pI3;
            newElementU.pointIndex[2] = pI1;
            newElementU.keepElement = true;

            //Info<< "after newElementU " << ' ' << eL.size() << endl;

            // allocation of the new element at the bottom part:
            elementInfo newElementB;
            newElementB.currentIndex = eL.size()+1;//eL[eL.size()-1].currentIndex + 1;
            newElementB.points[0] = newPoint.currentPoint;
            newElementB.points[1] = p1;
            newElementB.points[2] = p4;
            newElementB.pointIndex[0] = newPoint.currentIndex;
            newElementB.pointIndex[1] = pI1;
            newElementB.pointIndex[2] = pI4;
            newElementB.keepElement = true;


            //Info<< "after newElementB " << ' ' << eL.size() << endl;

            // transform currentElemet[I] and neighbElementI[c4III] to point to the newPoint 
            currentElemet.points[I] = newPoint.currentPoint;
            currentElemet.pointIndex[I] = newPoint.currentIndex;
            neighbElementI.points[c4III] = newPoint.currentPoint;
            neighbElementI.pointIndex[c4III] = newPoint.currentIndex;

            //Info<< "after transform currentElemet[I]  " << ' ' << eL.size() << endl;

            // update the neighboures of newElementU
            newElementU.elementIndex[0] = currentElemet.currentIndex;
            newElementU.elementIndex[1] = neighbElementIII.currentIndex;
            newElementU.elementIndex[2] = newElementB.currentIndex;

            //Info<< "after update the neighboures of newElementU " << ' ' << eL.size() << endl;

            // update the neighboures of newElementB
            newElementB.elementIndex[0] = newElementU.currentIndex;
            newElementB.elementIndex[1] = neighbElementILeft.currentIndex;
            newElementB.elementIndex[2] = neighbElementI.currentIndex;

            //Info<< "after update the neighboures of newElementB " << ' ' << eL.size() << endl;

            // update the neighboures of currentElemet and neighbElementI
            currentElemet.elementIndex[III] = newElementU.currentIndex;
            neighbElementI.elementIndex[c4III] = newElementB.currentIndex;

            //Info<< "after update the neighboures of currentElemet and neighbElementI " << ' ' << eL.size() << endl;

            // update the neighboures of neighbElementIII and neighbElementILeft
            neighbElementIII.elementIndex[c6II] = newElementU.currentIndex;
            neighbElementILeft.elementIndex[c7II] = newElementB.currentIndex;

            //Info<< "after update the neighboures of neighbElementIII and neighbElementILeft " << ' ' << eL.size() << endl;

            //updating the elOwner for p1
            //posInList = quickSelect(pL, pI1, 0, pL.size()-1);
            //pL[posInList].elOwner = neighbElementIII.currentIndex;
            pL[pI1].elOwner = neighbElementILeft.currentIndex;

            //Info << "--$$$---@@@@@28@@@@----neighbElementIII.elementIndex[c6I] " << neighbElementIII.elementIndex[0] << ' '
            //     << neighbElementIII.currentIndex << endl;
            //Info << "--$$$---@@@@@28@@@@----neighbElementILeft.elementIndex[c7I] " << neighbElementILeft.elementIndex[0] << ' '
            //     << neighbElementILeft.currentIndex << endl;

           // Info << "--$$$---@@@@@28@@@@----neighbElementIII.elementIndex[c6II] " << neighbElementIII.elementIndex[1] << ' '
           //      << neighbElementIII.currentIndex << endl;
           // Info << "--$$$---@@@@@28@@@@----neighbElementILeft.elementIndex[c7II] " << neighbElementILeft.elementIndex[1] << ' '
           //      << neighbElementILeft.currentIndex << endl;

           // Info << "--$$$---@@@@@28@@@@----neighbElementIII.elementIndex[c6III] " << neighbElementIII.elementIndex[2] << ' '
           //      << neighbElementIII.currentIndex << endl;
           // Info << "--$$$---@@@@@28@@@@----neighbElementILeft.elementIndex[c7III] " << neighbElementILeft.elementIndex[2] << ' '
           //      << neighbElementILeft.currentIndex << endl;

            pL.append(newPoint);
            eL.append(newElementU);
            eL.append(newElementB);

            Info<< "FFFFFFFFFFFffinal " << ' ' << elementI << endl;
       }//3
   }//2

   Info << "                 After refining " << numberOfELAdded << " element/elemens were added." << endl;
}//1
*/

/*
//-note:  func to refine the surface grid
template<class CloudType>
inline void Foam::FrontDynamic<CloudType>::frontRefining( DynamicList<pointData>& pL, DynamicList<elementInfo>& eL,
                            scalar minEdge, scalar maxEdge, scalar maxAspectRatio)
{//1
    label pI1, pI2, pI3, pI4;
    label c4I;
    //label c4II;
    label c4III;
    //label c5I
    //label c5II;
    //label c5III;
    //label c6I;
    label c6II;
    //label c6III;
    //label c7I;
    label c7II;
    //label c7III;
    //label c8I;
    //label c8II;
    //label c8III;
    label I, II, III;
    label f, cNI, cNII, cNIII;

    label numberOfELAdded = 0;
    
    label initialElListSize = eL.size();

    //forAll(eL, elementI)
    for(int elementI=0; elementI < initialElListSize; elementI++)      
    {//2
        elementInfo& currentElemet = eL[elementI];

            Info<< "AHey " << ' ' << elementI << endl;

        // order the lengths, s1 shortest, s3 longest 
        scalar s1 = mag( currentElemet.points[1] - currentElemet.points[0]);
        scalar s2 = mag( currentElemet.points[2] - currentElemet.points[1]);
        scalar s3 = mag( currentElemet.points[0] - currentElemet.points[2]);
        label n1 = 0;
        label n2 = 1;
        label n3 = 2;
            Info<< "s1/s2/s3 " << ' ' << s1 << ' ' << s2 << ' ' << s3 << endl;

        bool refinedElement = isElementRefined(currentElemet, minEdge, maxEdge, maxAspectRatio, s1, s2, s3, n1, n2, n3);

            Info<< "AAHey " << ' ' << refinedElement << currentElemet.keepElement << endl;

        if ( refinedElement && currentElemet.keepElement )
        {//3
            //Info<< "The element " << ' ' << currentElemet.currentIndex << ' ' << " will be refined." << endl;

            //pL.reserve(10*eL.size());
            //eL.reserve(10*eL.size());
                    //Info<< "bef rserving " << endl;
            //pL.reserve( pL.size()+1);
            //eL.reserve( eL.size()+2);

            numberOfELAdded = numberOfELAdded + 2;
             
            label largeEdge = n3;
                     
            // I, II, III are the point indexes of current element.
            // n1 is the index of shortest edge and p1.
            findPointsOrder( largeEdge, I, II, III);
            point p1 = currentElemet.points[I];//JJ
            point p2 = currentElemet.points[II];//JJ
            point p3 = currentElemet.points[III];//JJ
            pI1 = currentElemet.pointIndex[I];
            pI2 = currentElemet.pointIndex[II];
            pI3 = currentElemet.pointIndex[III];

            // find point 5, third point of neighbElementII
            //label posInList = quickSelect(eL, currentElemet.elementIndex[II], 0, eL.size()-1);
            //elementInfo& neighbElementII = eL[posInList];
            elementInfo& neighbElementII = eL[currentElemet.elementIndex[II]];
            for( f=0; f<=3; f++)
            {
                if( neighbElementII.pointIndex[f] != pI2 && neighbElementII.pointIndex[f] != pI3 )
                {
                    Info<< "f5 " << ' ' << f<< endl;
                    break;
                }
                if( f == 3)
                {
		    FatalErrorIn
		    (
			    " ^^^^^^^^^^^^^ we can not find such a point."
		    ) << abort(FatalError);
                }
            }            
            // c5I, c5II, c5III are the point indexes of neighbElementII.
            findPointsOrder( f, cNI, cNII, cNIII);
            point& p5 = neighbElementII.points[f];
            //c5I = cNI;
	    //c5II = cNII;
	    //c5III = cNIII;

            // find point 6, third point of neighbElementIII
            //posInList = quickSelect(eL, currentElemet.elementIndex[III], 0, eL.size()-1);
            //elementInfo& neighbElementIII = eL[posInList];
            elementInfo& neighbElementIII = eL[currentElemet.elementIndex[III]];
            for( f=0; f<=3; f++)
            {
                if( neighbElementIII.pointIndex[f] != pI1 && neighbElementIII.pointIndex[f] != pI3 )
                {
                    Info<< "f6 " << ' ' << f<< endl;
                    break;
                }
                if( f == 3)
                {
		    FatalErrorIn
		    (
			    " ^^^^^^^^^^^^^ we can not find such a point."
		    ) << abort(FatalError);
                }
            }
            findPointsOrder( f, cNI, cNII, cNIII);
            point& p6 = neighbElementIII.points[f];
	    //c6I = cNI;
	    c6II = cNII;
	    //c6III = cNIII;

            // find point 4, third point of neighbElementIII
            //posInList = quickSelect(eL, currentElemet.elementIndex[I], 0, eL.size()-1);
            //elementInfo& neighbElementI = eL[posInList];
            elementInfo& neighbElementI = eL[currentElemet.elementIndex[I]];
            for( f=0; f<=3; f++)
            {
                if( neighbElementI.pointIndex[f] != pI1 && neighbElementI.pointIndex[f] != pI2 )
                {
                    Info<< "f4 " << ' ' << f<< endl;
                    break;
                }
                if( f == 3)
                {
		    FatalErrorIn
		    (
			    " ^^^^^^^^^^^^^ we can not find such a point."
		    ) << abort(FatalError);
                }
            }
            findPointsOrder( f, cNI, cNII, cNIII);
            point& p4 = neighbElementI.points[f];
            pI4 = neighbElementI.pointIndex[f];
	    c4I = cNI;
	    //c4II = cNII;
	    c4III = cNIII;

            // find point 7, third point of neighbElementILeft
            //posInList = quickSelect(eL, neighbElementI.elementIndex[c4III], 0, eL.size()-1);
            //elementInfo& neighbElementILeft = eL[posInList];
            elementInfo& neighbElementILeft = eL[neighbElementI.elementIndex[c4III]];
            for( f=0; f<=3; f++)
            {
                if( neighbElementILeft.pointIndex[f] != pI1 && neighbElementILeft.pointIndex[f] != pI4 )
                {
                    Info<< "f7 " << ' ' << f<< endl;
                    break;
                }
                if( f == 3)
                {
		    FatalErrorIn
		    (
			    " ^^^^^^^^^^^^^ we can not find such a point."
		    ) << abort(FatalError);
                }
            }
            findPointsOrder( f, cNI, cNII, cNIII);
            point& p7 = neighbElementILeft.points[f];
	    //c7I = cNI;
	    c7II = cNII;
	    //c7III = cNIII;

            // find point 8, third point of neighbElementIRight
            //posInList = quickSelect(eL, neighbElementI.elementIndex[c4I], 0, eL.size()-1);
            //elementInfo& neighbElementIRight = eL[posInList];
            elementInfo& neighbElementIRight = eL[neighbElementI.elementIndex[c4I]];
            for( f=0; f<=3; f++)
            {
                if( neighbElementIRight.pointIndex[f] != pI2 && neighbElementIRight.pointIndex[f] != pI4 )
                {
                    Info<< "f8 " << ' ' << f<< endl;
                    break;
                }
                if( f == 3)
                {
		    FatalErrorIn
		    (
			    " ^^^^^^^^^^^^^ we can not find such a point."
		    ) << abort(FatalError);
                }
            }
            findPointsOrder( f, cNI, cNII, cNIII);
            point& p8 = neighbElementIRight.points[f];
	    //c8I = cNI;
	    //c8II = cNII;
	    //c8III = cNIII;

            // allocation of the new point:
            pointData newPoint;
            newPoint.currentIndex = pL.size();//pL[pL.size()-1].currentIndex + 1
            newPoint.currentPoint = 0.5 * ( p1 + p2) + 0.125 * ( p3 + p4) - 0.0625 * ( p5 + p6 + p7 + p8);
            newPoint.keepPoint = true;
            newPoint.elOwner = currentElemet.currentIndex;

            Info<< "   p1/p2/p3/p4/p5/p6/p7/p8: " << p1 << ' ' << p2 << ' ' << p3 << ' ' << p4 << ' '
                                                  << p5 << ' ' << p6 << ' ' << p7 << ' ' << p8 << endl;
            Info<< "newPoint for refining " << ' ' << newPoint.currentPoint << endl;

            // allocation of the new element at the upper part:
            elementInfo newElementU;
            newElementU.currentIndex = eL.size();//eL[eL.size()-1].currentIndex + 1;
            newElementU.points[0] = newPoint.currentPoint;
            newElementU.points[1] = p3;
            newElementU.points[2] = p1;
            newElementU.pointIndex[0] = newPoint.currentIndex;
            newElementU.pointIndex[1] = pI3;
            newElementU.pointIndex[2] = pI1;
            newElementU.keepElement = true;

            //Info<< "after newElementU " << ' ' << eL.size() << endl;

            // allocation of the new element at the bottom part:
            elementInfo newElementB;
            newElementB.currentIndex = eL.size()+1;//eL[eL.size()-1].currentIndex + 1;
            newElementB.points[0] = newPoint.currentPoint;
            newElementB.points[1] = p1;
            newElementB.points[2] = p4;
            newElementB.pointIndex[0] = newPoint.currentIndex;
            newElementB.pointIndex[1] = pI1;
            newElementB.pointIndex[2] = pI4;
            newElementB.keepElement = true;


            //Info<< "after newElementB " << ' ' << eL.size() << endl;

            // transform currentElemet[I] and neighbElementI[c4III] to point to the newPoint 
            currentElemet.points[I] = newPoint.currentPoint;
            currentElemet.pointIndex[I] = newPoint.currentIndex;
            neighbElementI.points[c4III] = newPoint.currentPoint;
            neighbElementI.pointIndex[c4III] = newPoint.currentIndex;

            //Info<< "after transform currentElemet[I]  " << ' ' << eL.size() << endl;

            // update the neighboures of newElementU
            newElementU.elementIndex[0] = currentElemet.currentIndex;
            newElementU.elementIndex[1] = neighbElementIII.currentIndex;
            newElementU.elementIndex[2] = newElementB.currentIndex;

            //Info<< "after update the neighboures of newElementU " << ' ' << eL.size() << endl;

            // update the neighboures of newElementB
            newElementB.elementIndex[0] = newElementU.currentIndex;
            newElementB.elementIndex[1] = neighbElementILeft.currentIndex;
            newElementB.elementIndex[2] = neighbElementI.currentIndex;

            //Info<< "after update the neighboures of newElementB " << ' ' << eL.size() << endl;

            // update the neighboures of currentElemet and neighbElementI
            currentElemet.elementIndex[III] = newElementU.currentIndex;
            neighbElementI.elementIndex[c4III] = newElementB.currentIndex;

            //Info<< "after update the neighboures of currentElemet and neighbElementI " << ' ' << eL.size() << endl;

            // update the neighboures of neighbElementIII and neighbElementILeft
            neighbElementIII.elementIndex[c6II] = newElementU.currentIndex;
            neighbElementILeft.elementIndex[c7II] = newElementB.currentIndex;

            //Info<< "after update the neighboures of neighbElementIII and neighbElementILeft " << ' ' << eL.size() << endl;

            //updating the elOwner for p1
            //posInList = quickSelect(pL, pI1, 0, pL.size()-1);
            //pL[posInList].elOwner = neighbElementIII.currentIndex;
            pL[pI1].elOwner = neighbElementILeft.currentIndex;

            Info << "--$$$---@@@@@28@@@@----neighbElementIII.elementIndex[c6I] " << neighbElementIII.elementIndex[0] << ' '
                 << neighbElementIII.currentIndex << endl;
            Info << "--$$$---@@@@@28@@@@----neighbElementILeft.elementIndex[c7I] " << neighbElementILeft.elementIndex[0] << ' '
                 << neighbElementILeft.currentIndex << endl;

            Info << "--$$$---@@@@@28@@@@----neighbElementIII.elementIndex[c6II] " << neighbElementIII.elementIndex[1] << ' '
                 << neighbElementIII.currentIndex << endl;
            Info << "--$$$---@@@@@28@@@@----neighbElementILeft.elementIndex[c7II] " << neighbElementILeft.elementIndex[1] << ' '
                 << neighbElementILeft.currentIndex << endl;

            Info << "--$$$---@@@@@28@@@@----neighbElementIII.elementIndex[c6III] " << neighbElementIII.elementIndex[2] << ' '
                 << neighbElementIII.currentIndex << endl;
            Info << "--$$$---@@@@@28@@@@----neighbElementILeft.elementIndex[c7III] " << neighbElementILeft.elementIndex[2] << ' '
                 << neighbElementILeft.currentIndex << endl;

            pL.append(newPoint);
            eL.append(newElementU);
            eL.append(newElementB);

            Info<< "FFFFFFFFFFFffinal " << ' ' << elementI << endl;
       }//3
   }//2

   Info << "                 After refining " << numberOfELAdded << " element/elemens were added." << endl;
}//1
*/

// note: func to removal front undulation based on TSUR method :
// Trapezoidal Sub-grid Undulations Removal. A front-tracking/front-capturing 
// method for the simulation of 3D multi-fluid flows with free surfaces, F.S. de Sousa et al,
// Journal of Computational Physics 198 (2004) 469–499
template<class CloudType>
inline void Foam::FrontDynamic<CloudType>::TSUR3D( label repeatNum, DynamicList<pointData>& pL, DynamicList<elementInfo>& eL)
{//1
    for( int repeat =0; repeat <= repeatNum-1; repeat++)  
    {//repeat
        // Vertex balance procedure
        forAll(pL, ptI) 
        {//2
       	    pointData& currentPt = pL[ptI];
            point& ver = currentPt.currentPoint;
	    DynamicList<point> cPoints = currentPt.connectedPoints;
                 
            // Normal vector calculation  
            // centre of the base for current point
            point p = vector::zero;
            forAll(cPoints, cPI)
            {
                p = p + cPoints[cPI];
            }
            p = p/cPoints.size();

            point xn = cPoints[cPoints.size()-1];
            point x0 = cPoints[0];
            vector normalSum = vector::zero;
            normalSum = normalSum + ( (xn - p) ^ (x0 - p) );
            for( int cPI = 0; cPI <= cPoints.size()-2; cPI++)
            {
                point xi = cPoints[cPI];
                point xii = cPoints[cPI+1];
                normalSum = normalSum + ( (xi - p) ^ (xii - p) );
            }
            vector n = normalSum/(mag(normalSum)+ROOTVSMALL);

            // calculation the V1 ,volume of the unitary polyhedron
            // (p+n, x0,..., xn, p)
            point pPlusN = p + n;
            scalar V1 = 0;
            vector Sf = 0.5 * (xn - p) ^ (x0 - p); 
            vector he = (p + xn + x0)/3.0 - pPlusN; 
            V1 = V1 + ( he & Sf)/3.0;
            for( int cPI = 0; cPI <= cPoints.size()-2; cPI++)
            {
                // volume of tetrahedron of p+n, p, xi, xii
                // he is a vector from p+n to the average of p, xi, and xii.
                point xi = cPoints[cPI];
                point xii = cPoints[cPI+1];

                Sf = 0.5 * (xi - p) ^ (xii - p); 
                he = (p + xi + xii)/3.0 - pPlusN; 
                V1 = V1 + ( he & Sf)/3.0;
            }
                  
            // calculation the V ,volume of the polyhedron
            // (ver, x0,..., xn, p)
            scalar V = 0;
            Sf = 0.5 * (xn - p) ^ (x0 - p); 
            he = (p + xn + x0)/3.0 - ver; 
            V = V + ( he & Sf)/3.0;
	    for( int cPI = 0; cPI <= cPoints.size()-2; cPI++)
            {
                // volume of tetrahedron of p+n, p, xi, xii
                // he is a vector from p+n to the average of p, xi, and xii.
                point xi = cPoints[cPI];
                point xii = cPoints[cPI+1];

                Sf = 0.5 * (xi - p) ^ (xii - p); 
                he = (p + xi + xii)/3.0 - ver; 
                V = V + ( he & Sf)/3.0;
             }

             scalar h = V/(V1+ROOTVSMALL);
             ver = p + h * n;
         }//2

         // Undulation removal procedure
         forAll(pL, ptI) 
         {//3
  	     pointData& curPt = pL[ptI];
	     DynamicList<point> cPoints = curPt.connectedPoints;
	     DynamicList<label> cPointsIndex = curPt.connectedPointsIndex;
			 
             point& ver1 = curPt.currentPoint;

	     // calculation of  p1 
  	     point p1 = vector::zero;
	     forAll(cPoints, cPI)
	     {
	         p1 = p1 + cPoints[cPI];
	     }
	     p1 = p1 /cPoints.size();

             forAll(cPoints, cPI)
             {//4
                 //posInList = quickSelect(pL, cPointsIndex[cPI], 0, pL.size()-1);
                 //pointData& otherPt = pL[posInList];
                 pointData& otherPt = pL[cPointsIndex[cPI]];
                 DynamicList<point> otherCPoints = otherPt.connectedPoints;

                 //-------------------------------------------------------------------
		 // calculation of V1(volume of the polyhedron), n1(normal vector) 

         	 // Normal vector calculation for curPt = ver1       
		 point x1n = cPoints[cPoints.size() -1];
		 point x10 = cPoints[0];
		 vector normalSum1 = vector::zero;
		 normalSum1 = normalSum1 + ( (x1n - p1) ^ (x10 - p1) );
		 for( int cPI = 0; cPI <= cPoints.size()-2; cPI++)
		 {
		     point x1i = cPoints[cPI];
		     point x1ii = cPoints[cPI+1];
		     normalSum1 = normalSum1 + ( (x1i - p1) ^ (x1ii - p1) );
		 }
		 vector n1 = normalSum1/(mag(normalSum1)+ROOTVSMALL); 

                 //----------------------------------------------------------------------
                 // calculation of V2(volume of the polyhedron), p2, n2(normal vector) 
                 point& ver2 = otherPt.currentPoint;

                 point p2 = vector::zero;
         	 forAll(otherCPoints, oCPI)
		 {
		     p2 = p2 + otherCPoints[oCPI];
                 }
		 p2 = p2 /otherCPoints.size();
                  
		 // Normal vector calculation for otherPt = ver2       
	       	 point x2n = otherCPoints[otherCPoints.size() -1];
	       	 point x20 = otherCPoints[0];
		 vector normalSum2 = vector::zero;
		 normalSum2 = normalSum2 + ( (x2n - p2) ^ (x20 - p2) );
		 for( int oCPI = 0; oCPI <= otherCPoints.size()-2; oCPI++)
		 {
		     point x2i = otherCPoints[oCPI];
		     point x2ii = otherCPoints[oCPI+1];
		     normalSum2 = normalSum2 + ( (x2i - p2) ^ (x2ii - p2) );
		 }
		 vector n2 = normalSum2/(mag(normalSum2)+ROOTVSMALL); 

                 //----------------------------------------------------------------------
                 //----------------------------------------------------------------------
		 vector n = (n1 + n2)/mag(n1 + n2);
                 
		 vector m = 0.5 * (ver1 + ver2);
                 
		 scalar h1 = mag( (ver1 - m) & n );
		 scalar h2 = mag( (ver2 - m) & n );
                 
                 // new definitions for p1 and p2
		 p1 = ver1 - (h1 * n);
		 p2 = ver2 - (h2 * n);

                 //----------------------------------------------------------------------
 		 // calculation the V1 ,volume of the polyhedron
		 // ver1, x10,..., x1n, p1)
		 scalar V1 = 0;
		 vector Sf1 = 0.5 * (x1n - p1) ^ (x10 - p1); 
		 vector he1 = (p1 + x1n + x10)/3.0 - ver1; 
		 V1 = V1 + mag( he1 & Sf1)/3.0;
         	 for( int cPI = 0; cPI <= cPoints.size()-2; cPI++)
		 {
		     point x1i = cPoints[cPI];
		     point x1ii = cPoints[cPI+1];

		     //volume of tetrahedron of ver1, p1, x1i, x1ii
		     //he1 is a vector from ver1 to the average of p1, x1i, and x1ii.
		     Sf1 = 0.5 * (x1i - p1) ^ (x1ii - p1); 
		     he1 = (p1 + x1i + x1ii)/3.0 - ver1; 

		     V1 = V1 + mag( he1 & Sf1)/3.0;
		 }

                 //----------------------------------------------------------------------
	         // calculation the V2 ,volume of the polyhedron
	         // (ver, x20,..., x2n, p2)
	         scalar V2 = 0;
                 vector Sf2 = 0.5 * (x2n - p2) ^ (x20 - p2); 
                 vector he2 = (p2 + x2n + x20)/3.0 - ver2; 
                 V2 = V2 + mag( he2 & Sf2)/3.0;
    		 for( int oCPI = 0; oCPI <= otherCPoints.size()-2; oCPI++)
	         {
	             point x2i = otherCPoints[oCPI];
	             point x2ii = otherCPoints[oCPI+1];

		     // volume of tetrahedron of ver2, p2, xi, xii
		     // h is a vector from p+n to the average of p, xi, and xii.
                     vector Sf2 = 0.5 * (x2i - p2) ^ (x2ii - p2); 
	             vector he2 = (p2 + x2i + x2ii)/3.0 - ver2; 
		      
                     V2 = V2 + mag( he2 & Sf2)/3.0;
		 }

                 //----------------------------------------------------------------------
	         scalar h = (V1 + V2) / ( mag( V1/(h1 + ROOTVSMALL) ) + mag( V2/(h2 + ROOTVSMALL) ) );
	         ver1 = p1 + h * n;
 	         ver2 = p2 + h * n;
             }//4
         }//3

    }//repeat
   
    pointToElementMapping(pL, eL);
}//1

// note: func to smooth the front mesh, algorithm III.
// Volume Conserving Smoothing for Piecewise Linear Curves, Surfaces, and Triple Lines 
//A. Kuprat et al, Journal of Computational Physics 172, 99–118 (2001)
template<class CloudType>
inline void Foam::FrontDynamic<CloudType>::volumeConservingSmoothingIII( label repeatNum, DynamicList<pointData>& pL, DynamicList<elementInfo>& eL)
{//1

    for( int repeat =0; repeat <= repeatNum-1; repeat++)  
    {//repeat

        // Vertex balance procedure
        forAll(pL, ptI) 
        {//2
       	    pointData& currentPt = pL[ptI];
            point& x = currentPt.currentPoint;
	    DynamicList<point>& cPoints = currentPt.connectedPoints;
                 
            // calculation of the 
            point xs = vector::zero;
            forAll(cPoints, cPI)
            {
                xs = xs + cPoints[cPI];
            }
            xs = xs/cPoints.size();

            point xn = cPoints[cPoints.size()-1];
            point x0 = cPoints[0];
            vector normalSum = vector::zero;
            normalSum = normalSum + ( (xn - x) ^ (x0 - x) );
            for( int cPI = 0; cPI <= cPoints.size()-2; cPI++)
            {
                point xi = cPoints[cPI];
                point xii = cPoints[cPI+1];
                normalSum = normalSum + ( (xi - x) ^ (xii - x) );
            }
            vector n = normalSum/(mag(normalSum)+ROOTVSMALL);

            // calculation of the dxs
            vector dxs = xs - x;

            // calculation of the smoothed x
            x = x + dxs - ( dxs & n ) * n;
         }//2
    }//repeat
    pointToElementMapping(pL, eL);
}//1

// note: func to smooth the front mesh, algorithm IV.
// Volume Conserving Smoothing for Piecewise Linear Curves, Surfaces, and Triple Lines 
//A. Kuprat et al, Journal of Computational Physics 172, 99–118 (2001)
template<class CloudType>
inline void Foam::FrontDynamic<CloudType>::volumeConservingSmoothingIV( label repeatNum, DynamicList<pointData>& pL, DynamicList<elementInfo>& eL)
{//1
    //scalar volDiff = (calcBubbleVolume(eL) - initialBubbleVolume_)/eL.size();
    scalar volDiff = 0.0; //(calcBubbleVolume(eL) - initialBubbleVolume_)/pL.size();

    scalar omega = 1.0;
    label f;

    for( int repeat =0; repeat <= repeatNum-1; repeat++)  
    {//repeat
         // Undulation removal procedure
         forAll(pL, ptI) 
         {//2
  	     pointData& curPt = pL[ptI];
	     DynamicList<point> cPoints = curPt.connectedPoints;
	     DynamicList<label> cPointsIndex = curPt.connectedPointsIndex;
             point& x2 = curPt.currentPoint;

             // calculation of the A2
             point xn = cPoints[cPoints.size()-1];
             point x0 = cPoints[0];
             vector A2 = vector::zero;
             A2 = A2 + ( (xn - x2) ^ (x0 - x2) );
             for( int cPI = 0; cPI <= cPoints.size()-2; cPI++)
             {
                 point xi = cPoints[cPI];
                 point xii = cPoints[cPI+1];
                 A2 = A2 + ( (xi - x2) ^ (xii - x2) );
             }

             forAll(cPoints, cPI)
             {//4
                 pointData& otherPt = pL[cPointsIndex[cPI]];
                 DynamicList<point> otherCPoints = otherPt.connectedPoints;
                 point& x1 = otherPt.currentPoint;

                 // calculation of the A1
                 xn = otherCPoints[otherCPoints.size()-1];
                 x0 = otherCPoints[0];
                 vector A1 = vector::zero;
                 A1 = A1 + ( (xn - x1) ^ (x0 - x1) );
                 for( int oCPI = 0; oCPI <= otherCPoints.size()-2; oCPI++)
                 {
                     point xi = otherCPoints[oCPI];
                     point xii = otherCPoints[oCPI+1];
                     A1 = A1 + ( (xi - x1) ^ (xii - x1) );
                 }

                 // finding e22 for v
                 label up, down;
                 for( f=0; f <= cPoints.size()-1; f++)
                 {
                     if( curPt.connectedPointsIndex[f] == otherPt.currentIndex )
                     {
                         break;
                     }
                 }
                 if ( f == 0)
                 {
                     up = 1;
                     down = cPoints.size()-1; 
                 }
                 else if ( f == (cPoints.size()-1) )
                 {
                     up = 0;
                     down = cPoints.size()-2; 
                 }
                 else
                 {
                     up = f+1;
                     down = f-1; 
                 }
                 vector e2n2 = cPoints[down] - x2;
                 vector e22 = cPoints[up] - x2;
                 vector v = e2n2 - e22;

                 // calculation of the dx1s and dx2s
                 vector dx1s = vector::zero;
                 vector dx2s = vector::zero;

                 if( undulationRemovalingType_ == "volumeConservingSmoothingIV" )
                 {
                     // calculation of the sumX2
                     vector sumX2 = vector::zero;
                     for( int oCPI = 0; oCPI <= cPoints.size()-1; oCPI++)
                     {
                         sumX2 = sumX2 + cPoints[oCPI];
                     }
                     sumX2 = sumX2 - x1; 

                     // calculation of the sumX1
                     vector sumX1 = vector::zero;
                     for( int oCPI = 0; oCPI <= otherCPoints.size()-1; oCPI++)
                     {
                         sumX1 = sumX1 + otherCPoints[oCPI];
                     }
                     sumX1 = sumX1 - x2;
 
                     // calculation of the x1s and x2s
                     scalar n1 = otherCPoints.size();
                     scalar n2 = cPoints.size();
                     vector x1s = 1./( n1*n2-1.) * (sumX2 + n2 * sumX1);
                     vector x2s = (1./n2) * (x1s + sumX2);

                     // calculation of the dx1s and dx2s
                     dx1s = omega * (x1s - x1);
                     dx2s = omega * (x2s - x2);
                 }
                 if( undulationRemovalingType_ == "combined" )
                 {
                     // calculation of the dx1s and dx2s
                     dx1s = omega * (smoothTo(otherPt) - x1);
                     dx2s = omega * (smoothTo(curPt) - x2);
                 }                

                 // calculation of the A
                 vector A = A1 + A2 + (v ^ (dx1s - dx2s));

                 // calculation of the normal,n
                 vector n = A/(mag(A)+ROOTVSMALL);

                 // calculation of the h
                 scalar h = -( (dx1s & A1) + (dx2s & A2) + ( dx2s & (v ^ dx1s) ) )/( (n & A) + volDiff + ROOTVSMALL);

                 // calculation of the new x2 and x1
                 x2 = x2 + dx2s + h * n;
                 x1 = x1 + dx1s + h * n;
             }//4
         }//2
    }//repeat
    pointToElementMapping(pL, eL);
}//1

//- note: func to smooth a pointData var.
template<class CloudType>
inline Foam::vector Foam::FrontDynamic<CloudType>::smoothTo(pointData currentPt)
{
            point ver = currentPt.currentPoint;
	    DynamicList<point> cPoints = currentPt.connectedPoints;
                 
            // Normal vector calculation  
            // centre of the base for current point
            point p = vector::zero;
            forAll(cPoints, cPI)
            {
                p = p + cPoints[cPI];
            }
            p = p/cPoints.size();

            point xn = cPoints[cPoints.size()-1];
            point x0 = cPoints[0];
            vector normalSum = vector::zero;
            normalSum = normalSum + ( (xn - p) ^ (x0 - p) );
            for( int cPI = 0; cPI <= cPoints.size()-2; cPI++)
            {
                point xi = cPoints[cPI];
                point xii = cPoints[cPI+1];
                normalSum = normalSum + ( (xi - p) ^ (xii - p) );
            }
            vector n = normalSum/(mag(normalSum)+ROOTVSMALL);

            // calculation the V1 ,volume of the unitary polyhedron
            // (p+n, x0,..., xn, p)
            point pPlusN = p + n;
            scalar V1 = 0;
            vector Sf = 0.5 * (xn - p) ^ (x0 - p); 
            vector he = (p + xn + x0)/3.0 - pPlusN; 
            V1 = V1 + ( he & Sf)/3.0;
            for( int cPI = 0; cPI <= cPoints.size()-2; cPI++)
            {
                // volume of tetrahedron of p+n, p, xi, xii
                // he is a vector from p+n to the average of p, xi, and xii.
                point xi = cPoints[cPI];
                point xii = cPoints[cPI+1];

                Sf = 0.5 * (xi - p) ^ (xii - p); 
                he = (p + xi + xii)/3.0 - pPlusN; 
                V1 = V1 + ( he & Sf)/3.0;
            }
                  
            // calculation the V ,volume of the polyhedron
            // (ver, x0,..., xn, p)
            scalar V = 0;
            Sf = 0.5 * (xn - p) ^ (x0 - p); 
            he = (p + xn + x0)/3.0 - ver; 
            V = V + ( he & Sf)/3.0;
	    for( int cPI = 0; cPI <= cPoints.size()-2; cPI++)
            {
                // volume of tetrahedron of p+n, p, xi, xii
                // he is a vector from p+n to the average of p, xi, and xii.
                point xi = cPoints[cPI];
                point xii = cPoints[cPI+1];

                Sf = 0.5 * (xi - p) ^ (xii - p); 
                he = (p + xi + xii)/3.0 - ver; 
                V = V + ( he & Sf)/3.0;
             }

             scalar h = V/(V1+ROOTVSMALL);
             ver = p + h * n;
             return ver;
}

//- note: func to calculate the volume enclosed by a front.
template<class CloudType>
inline Foam::scalar Foam::FrontDynamic<CloudType>::calcBubbleVolume( DynamicList<elementInfo> eL)
{
    scalar bubbleVolume = 0;
    scalar curTetrahedralVolume = 0;
    vector elementCenter = vector::zero;

    forAll(eL, elementI)  
    {
        elementInfo& elD = eL[elementI];

        vector e01 = ( elD.points[1] - elD.points[0]); 
        vector e02 = ( elD.points[2] - elD.points[0]);
        vector Sf = 0.5 * (e01 ^ e02);

        elementCenter = ( elD.points[0] + elD.points[1] + elD.points[2])/3.0;

        curTetrahedralVolume = ( elementCenter & Sf )/3.0;

        bubbleVolume = bubbleVolume + curTetrahedralVolume;
    }

    return mag(bubbleVolume);
}


//- note: func to calculate the bubble surface or front surface.
template<class CloudType>
inline Foam::scalar Foam::FrontDynamic<CloudType>::calcBubbleSurface( DynamicList<elementInfo> eL)
{
    scalar bubbleSurface = 0;
    forAll(eL, elementI)  
    {
        elementInfo& elD = eL[elementI];
        vector e01 = ( elD.points[1] - elD.points[0]); 
        vector e02 = ( elD.points[2] - elD.points[0]);
        vector Sf = 0.5 * (e01 ^ e02);
        bubbleSurface = bubbleSurface + mag(Sf);
    }
    return bubbleSurface;
}

//- note: func to calculate the Peskin distribution for structured Mesh
template<class CloudType>
inline Foam::scalar Foam::FrontDynamic<CloudType>::PeskinD(scalar r)
{
   if( mag(r) <= 2 )
   {
      return (1.0/4.0)*( 1.0 + cos( mathematical::piByTwo * r) ); 
   }
   else
   {
      return 0.0; 
   }
}

//- note: func to calculate the window function for structured Mesh 
template<class CloudType>
inline Foam::scalar Foam::FrontDynamic<CloudType>::windowFunc(vector del, scalar h)
{
    scalar dX = del.x(); 
    scalar dY = del.y();
    scalar dZ = del.z();
    return PeskinD(dX) * PeskinD(dY) * PeskinD(dZ); 
}

//- note: func to calculate the window function for unstructured Mesh
template<class CloudType>
inline Foam::scalar Foam::FrontDynamic<CloudType>::windowFuncBar(vector del, scalar h, List<scalar> beta)
{
    scalar dX = del.x();
    scalar dY = del.y();
    scalar dZ = del.z();
    return windowFunc(del, h) * ( beta[0] + beta[1] * dX + beta[2] * dY + beta[3] * dZ );
}

//- note: func to calculate the surface tension
template<class CloudType>
inline Foam::scalar Foam::FrontDynamic<CloudType>::mPQR(vector del, scalar h, label p, label q, label r)
{
    scalar dX = del.x();
    scalar dY = del.y();
    scalar dZ = del.z();
    return pow(dX,p) * pow(dY,q) * pow(dZ,r) *  windowFunc(del, h);
}

//- note: func to calculate the element surface area formed
//  by points connected to the current front point
template<class CloudType>
inline Foam::vector Foam::FrontDynamic<CloudType>::calcTheSurfaceAreaAtpoint( pointData ptD)
{
    vector surfaceArea = vector::zero; 
    DynamicList<point> cPoints = ptD.connectedPoints;

    vector x1 = ptD.currentPoint;//currentPoint
    vector x2 = cPoints[0];//pt1
    vector x3 = cPoints[1];//pt2

    vector e12 = ( x2 - x1); 
    vector e13 = ( x3 - x1);
    vector areaEnclosedByX1X2X3 = 0.5 * (e12 ^ e13); 
    surfaceArea =  surfaceArea + areaEnclosedByX1X2X3;
 
    for( int cPI = 1; cPI <= cPoints.size()-2; cPI++)
    {
        x2 = cPoints[cPI];//pt1
        x3 = cPoints[cPI + 1];//pt2

        e12 = ( x2 - x1); 
        e13 = ( x3 - x1);
        areaEnclosedByX1X2X3 = 0.5 * (e12 ^ e13); 
        surfaceArea =  surfaceArea + areaEnclosedByX1X2X3;
    }
    return surfaceArea/cPoints.size();
}

//- note: func to calculate the surface tension force at each point
template<class CloudType>
inline Foam::vector Foam::FrontDynamic<CloudType>::calcSurfaceTensionForceAtpoint( pointData ptD)
{
    vector surfaceTensionForce = vector::zero; 
    DynamicList<point> cPoints = ptD.connectedPoints;

    vector x1 = ptD.currentPoint;//currentPoint
    vector x3 = cPoints[0];//pt1
    vector x2 = cPoints[1];//pt2

    vector delX12 = (x1 - x2); 
    vector delX13 = (x1 - x3); 
    vector delX32 = (x3 - x2); 

    vector ne = (delX13 ^ delX12)/mag(delX13 ^ delX12); // this is corrected. because in Trigvasson book,
                                                        // the normal direction is calculated incorrectly.

    surfaceTensionForce =  surfaceTensionForce + 0.5 * (ne ^ delX32);
 
    for( int cPI = 1; cPI <= cPoints.size()-2; cPI++)
    {
        x3 = cPoints[cPI];//pt1
        x2 = cPoints[cPI + 1];//pt2

        delX12 = (x1 - x2); 
        delX13 = (x1 - x3); 
        delX32 = (x3 - x2); 
        ne = (delX13 ^ delX12)/mag(delX13 ^ delX12); // this is corrected. because in Trigvasson book,
                                                     // the normal direction is calculated incorrectly. 
        surfaceTensionForce =  surfaceTensionForce + 0.5 * (ne ^ delX32);
    }
    return surfaceTensionForce;
}

//- note: func to calculate the surface area of each face
template<class CloudType>
inline Foam::vector Foam::FrontDynamic<CloudType>::calcElementSurfaceArea( elementInfo elD)
{
    vector e10 = ( elD.points[1] - elD.points[0]); 
    vector e20 = ( elD.points[2] - elD.points[0]);
    return 0.5 * (e10 ^ e20);
}

//note: func to print the front data for all times.
template<class CloudType>
inline void Foam::FrontDynamic<CloudType>::printFrontsForAllTimes()
{
    fileName outfileName = "printFrontsForAllTimesPLT";

    // This is the directory keeping the bubble data. 
    fileName FTResult = "FTResult";

    fileName outfileNamePath=this->owner().db().time().rootPath()/this->owner().db().time().globalCaseName()/FTResult/outfileName;

    ofstream outFile(outfileNamePath.c_str(), ios::app);
    ostream& os = outFile;

    label bubbleNum = bDataL_.size();

    const scalar time = this->owner().db().time().value();

    // preparing the zone name according to the time value.
    OStringStream oss;
    oss << "BubbleAtTime=" << time;
    fileName zoneName = oss.str();

    label allElCounter = 0;
    label allPtCounter = 0;
    for(int bI = 0; bI < bubbleNum ; bI++)
    {//0
        bubbleData bData;
        DynamicList<pointData>& pL = bDataL_[bI].pL;
        DynamicList<elementInfo>& eL = bDataL_[bI].eL;
        allElCounter = allElCounter + eL.size();
        allPtCounter = allPtCounter + pL.size();
    }//0

    // Write header
    os  << " VARIABLES = \"X\" \"Y\" \"Z\",\"NUM\" "
	<< "\n ZONE T=\"" << zoneName << "\""
	<< "\n N=        " << allPtCounter
	<< " E=        " << allElCounter 
	<< " ZONETYPE=FETriangle" 
	<< "\n DATAPACKING=POINT \n";
      
    for(int bI = 0; bI < bubbleNum ; bI++)
    {//0
        DynamicList<pointData> pL = bDataL_[bI].pL;
        // Write vertex coords
        forAll(pL, pI)
        {
            vector curPoint = pL[pI].currentPoint;
	    os  << curPoint.x() << token::SPACE
	        << curPoint.y() << token::SPACE
	        << curPoint.z() << token::SPACE << bI << " \n";
        }

    }//0

    for(int bI = 0; bI < bubbleNum ; bI++)
    {//0
        DynamicList<elementInfo> eL = bDataL_[bI].eL;
        forAll(eL, eI)
        {
	    os  << eL[eI].pointIndex[0]+1 << token::SPACE               
	        << eL[eI].pointIndex[1]+1 << token::SPACE
	        << eL[eI].pointIndex[2]+1 << token::SPACE << " \n";
        }
    }//0

    os  << " \n";

    outFile.close();
}

//note: func to write data of all bubbles.
template<class CloudType>
inline void Foam::FrontDynamic<CloudType>::writeDataOfTheFrontsAtTheOutputTime()
{
    fileName outfileName = "dataOfTheFrontsAtTheCurrentTime";

    // This is the directory keeping the bubble data. 
    fileName FTResult = "FTResult";

    fileName outfileNamePath=this->owner().db().time().rootPath()/this->owner().db().time().globalCaseName()/FTResult/outfileName;

    ofstream outFile(outfileNamePath.c_str(), ios::out);
    ostream& os = outFile;

    label bubbleNum = bDataL_.size();

    const scalar time = this->owner().db().time().value();

    os  << "# Time = " << time << " \n";
    os  << "#" << " \n";
    os  << "#" << " \n";
    os  << "# The number of bubbles is as follows: " << " \n";
    os  << "  " << bubbleNum << " \n";

    for(int bI = 0; bI < bubbleNum ; bI++)
    {//0
        DynamicList<pointData> pL = bDataL_[bI].pL;
        DynamicList<elementInfo> eL =  bDataL_[bI].eL;

        os  << "# Data for bubble number: " << bI << " \n";
        os  << "# The number of points and elements are as follows: " << " \n";
        os  << "  " << pL.size() << ' ' << eL.size() << " \n";
        os  << "# The outerFluidDensity/density/outerFluidViscosity/viscosity/surfaceTensionCoeff/Eotvos/Morton are as follows: " << " \n";
        os  << "  " << bDataL_[bI].outerFluidDensity << ' ' << bDataL_[bI].density << ' '
                    << bDataL_[bI].outerFluidViscosity << ' ' << bDataL_[bI].viscosity << ' '
                    << bDataL_[bI].surfaceTensionCoeff << ' ' << bDataL_[bI].EotvosNumber  << ' ' << bDataL_[bI].MortonNumber << " \n";
        os  << "# The initial bubble diameter is as follows: " << " \n";
        os  << "  " << bDataL_[bI].sphereBubbleDiameter << " \n";
        os  << "# X, Y, and Z components of point positions are as follows: " << " \n";
        forAll(pL, pI)
        {
            vector curPoint = pL[pI].currentPoint;
	    os  << curPoint.x() << token::SPACE
	        << curPoint.y() << token::SPACE
	        << curPoint.z() << token::SPACE << " \n";
        }
        os  << "# all point indexes of each element are as follows: " << "\n";
        forAll(eL, eI)
        {
	    os  << eL[eI].pointIndex[0] << token::SPACE               
	        << eL[eI].pointIndex[1] << token::SPACE
	        << eL[eI].pointIndex[2] << token::SPACE << " \n";
        }

    }//0

    outFile.close();
}

//note: func to write data of all bubbles  .
template<class CloudType>
inline void Foam::FrontDynamic<CloudType>::writeDataOfTheFrontsAtTheOutputTimeIntimeName()
{

    fileName outfileName = "dataOfTheFrontsAtTheCurrentTime";

    mkDir
    ( 
        this->owner().db().time().rootPath()/this->owner().db().time().globalCaseName()/this->owner().db().time().timeName()
    ); 

    // This is the directory keeping the bubble data. 
    fileName outfileNamePath = this->owner().db().time().rootPath()/this->owner().db().time().globalCaseName()/this->owner().db().time().timeName()/outfileName;

    ofstream outFile(outfileNamePath.c_str(), ios::out);
    ostream& os = outFile;

    label bubbleNum = bDataL_.size();

    const scalar time = this->owner().db().time().value();

    os  << "# Time = " << time << " \n";
    os  << "#" << " \n";
    os  << "#" << " \n";
    os  << "# The number of bubbles is as follows: " << " \n";
    os  << "  " << bubbleNum << " \n";

    for(int bI = 0; bI < bubbleNum ; bI++)
    {//0
        DynamicList<pointData> pL = bDataL_[bI].pL;
        DynamicList<elementInfo> eL =  bDataL_[bI].eL;

        os  << "# Data for bubble number: " << bI << " \n";
        os  << "# The number of points and elements are as follows: " << " \n";
        os  << "  " << pL.size() << ' ' << eL.size() << " \n";
        os  << "# The outerFluidDensity/density/outerFluidViscosity/viscosity/surfaceTensionCoeff/Eotvos/Morton are as follows: " << " \n";
        os  << "  " << bDataL_[bI].outerFluidDensity << ' ' << bDataL_[bI].density << ' '
                    << bDataL_[bI].outerFluidViscosity << ' ' << bDataL_[bI].viscosity << ' '
                    << bDataL_[bI].surfaceTensionCoeff << ' ' << bDataL_[bI].EotvosNumber  << ' ' << bDataL_[bI].MortonNumber << " \n";
        os  << "# The initial bubble diameter is as follows: " << " \n";
        os  << "  " << bDataL_[bI].sphereBubbleDiameter << " \n";
        os  << "# X, Y, and Z components of point positions are as follows: " << " \n";
        forAll(pL, pI)
        {
            vector curPoint = pL[pI].currentPoint;
	    os  << curPoint.x() << token::SPACE
	        << curPoint.y() << token::SPACE
	        << curPoint.z() << token::SPACE << " \n";
        }
        os  << "# all point indexes of each element are as follows: " << "\n";
        forAll(eL, eI)
        {
	    os  << eL[eI].pointIndex[0] << token::SPACE               
	        << eL[eI].pointIndex[1] << token::SPACE
	        << eL[eI].pointIndex[2] << token::SPACE << " \n";
        }

    }//0

    outFile.close();
}

//note: func to write data of all bubbles.
template<class CloudType>
inline void Foam::FrontDynamic<CloudType>::readDataOfTheFrontsAtTheCurrentTime()
{
    Info << "\n---------> At the time: " << this->owner().db().time().value()
         <<" , Now all data of the front is reading to continue solving ....... " << endl;

    fileName infileName = "dataOfTheFrontsAtTheCurrentTime";

    // This is the directory keeping the bubble data. 
    fileName FTResult = "FTResult";

    fileName InfileNamePath=this->owner().db().time().rootPath()/this->owner().db().time().globalCaseName()/FTResult/infileName;

    IFstream inFile(InfileNamePath);

    label bubbleNum;
    pointData ptI; //tempoaray pointData var
    elementInfo elI; //tempoaray elementInfo var
    label ptNum, elNum; //tempoaray vars
    scalar x, y, z; //tempoaray vars
    label poI, poII, poIII; //tempoaray vars

    string line = getLineNoComment(inFile);
    IStringStream lineStream(line);
    lineStream >> bubbleNum;

    List<fileName> printAllbubbles(bubbleNum);
    forAll(printAllbubbles,pAB)
    { 
        // note, this way is for intrinsic C++ utility
	//std::ostringstream oss;
	//oss << "b" << pAB << ".plt";
	//printAllbubbles[pAB] = oss.str();

        // note, this way is special for OpenFOAM using the OStringStream object
        OStringStream oss;
	oss << "bubble" << pAB << "_AtCurrentTime" << ".plt";
        fileName fName = oss.str();
	printAllbubbles[pAB] =  this->owner().db().time().rootPath()/this->owner().db().time().globalCaseName()/FTResult/fName;
    }

    for(int bI = 0; bI < bubbleNum ; bI++)
    {//00
        bubbleData bData;
        DynamicList<pointData>& pL = bData.pL;
        DynamicList<elementInfo>& eL = bData.eL;

        pL.reserve( 100000);//JJ
        eL.reserve( 200000);//JJ
        
        string line = getLineNoComment(inFile);
        IStringStream lineStream(line);
        lineStream >> ptNum >> elNum;

        string line2 = getLineNoComment(inFile);
        IStringStream lineStream2(line2);
        lineStream2 >> bData.outerFluidDensity >> bData.density
                    >> bData.outerFluidViscosity >> bData.viscosity
                    >> bData.surfaceTensionCoeff >> bData.EotvosNumber >> bData.MortonNumber;

        string line3 = getLineNoComment(inFile);
        IStringStream lineStream3(line3);
        lineStream3 >> bData.sphereBubbleDiameter;

        for(int I = 0; I < ptNum ; I++)
        {//11
            string line = getLineNoComment(inFile);
            {
                IStringStream lineStream(line);
                lineStream >> x >> y >> z;
            }    
            ptI.currentPoint = point(x, y, z);
            ptI.currentIndex = I;
            ptI.keepPoint = true;
            pL.append(ptI);
        }//11

        for(int II = 0; II < elNum ; II++)
        {//12
            string line = getLineNoComment(inFile);
            {
                IStringStream lineStream(line);
                lineStream >> poI >> poII >> poIII;
            }   
            elI.keepElement = true;
            elI.currentIndex = II;
            elI.pointIndex[0] = poI;
            elI.pointIndex[1] = poII;
            elI.pointIndex[2] = poIII;
            eL.append(elI);
        }//13

        startPtIndexFromZero(pL,eL);
        findNighboures(pL,eL);
	pointToElementMapping(pL,eL);
	ptConnectedPoints(pL,eL);
        printInitialFronts(printAllbubbles[bI], pL, eL);
        bDataL_.append(bData);
    }//00
    initialAllocatingTheCPTList();
}

//note: func to write data of all bubbles.
template<class CloudType>
inline void Foam::FrontDynamic<CloudType>::readDataOfTheFrontsAtTheCurrentTimeIntimeName()
{
    Info << "\n---------> At the time: " << this->owner().db().time().value()
         <<" , Now all data of the front is reading to continue solving ....... " << endl;

    fileName infileName = "dataOfTheFrontsAtTheCurrentTime";

    // This is the directory keeping the bubble data.
    fileName FTResult = "FTResult"; 
    fileName InfileNamePath=this->owner().db().time().rootPath()/this->owner().db().time().globalCaseName()/this->owner().db().time().timeName()/infileName;

    IFstream inFile(InfileNamePath);

    label bubbleNum;
    pointData ptI; //tempoaray pointData var
    elementInfo elI; //tempoaray elementInfo var
    label ptNum, elNum; //tempoaray vars
    scalar x, y, z; //tempoaray vars
    label poI, poII, poIII; //tempoaray vars

    string line = getLineNoComment(inFile);
    IStringStream lineStream(line);
    lineStream >> bubbleNum;

    List<fileName> printAllbubbles(bubbleNum);
    forAll(printAllbubbles,pAB)
    { 
        // note, this way is for intrinsic C++ utility
	//std::ostringstream oss;
	//oss << "b" << pAB << ".plt";
	//printAllbubbles[pAB] = oss.str();

        // note, this way is special for OpenFOAM using the OStringStream object
        OStringStream oss;
	oss << "bubble" << pAB << "_AtCurrentTime" << ".plt";
        fileName fName = oss.str();
	printAllbubbles[pAB] =  this->owner().db().time().rootPath()/this->owner().db().time().globalCaseName()/FTResult/fName;
    }

    for(int bI = 0; bI < bubbleNum ; bI++)
    {//00
        bubbleData bData;
        DynamicList<pointData>& pL = bData.pL;
        DynamicList<elementInfo>& eL = bData.eL;

        pL.reserve( 100000);//JJ
        eL.reserve( 200000);//JJ
        
        string line = getLineNoComment(inFile);
        IStringStream lineStream(line);
        lineStream >> ptNum >> elNum;

        string line2 = getLineNoComment(inFile);
        IStringStream lineStream2(line2);
        lineStream2 >> bData.outerFluidDensity >> bData.density
                    >> bData.outerFluidViscosity >> bData.viscosity
                    >> bData.surfaceTensionCoeff >> bData.EotvosNumber >> bData.MortonNumber;

        string line3 = getLineNoComment(inFile);
        IStringStream lineStream3(line3);
        lineStream3 >> bData.sphereBubbleDiameter;

        for(int I = 0; I < ptNum ; I++)
        {//11
            string line = getLineNoComment(inFile);
            {
                IStringStream lineStream(line);
                lineStream >> x >> y >> z;
            }    
            ptI.currentPoint = point(x, y, z);
            ptI.currentIndex = I;
            ptI.keepPoint = true;
            pL.append(ptI);
        }//11

        for(int II = 0; II < elNum ; II++)
        {//12
            string line = getLineNoComment(inFile);
            {
                IStringStream lineStream(line);
                lineStream >> poI >> poII >> poIII;
            }   
            elI.keepElement = true;
            elI.currentIndex = II;
            elI.pointIndex[0] = poI;
            elI.pointIndex[1] = poII;
            elI.pointIndex[2] = poIII;
            eL.append(elI);
        }//13

        startPtIndexFromZero(pL,eL);
        findNighboures(pL,eL);
	pointToElementMapping(pL,eL);
	ptConnectedPoints(pL,eL);
        printInitialFronts(printAllbubbles[bI], pL, eL);
        bDataL_.append(bData);
    }//00
    initialAllocatingTheCPTList();
    allCellsInEeachMasket_.setSize(bDataL_.size());
  //allCellsNearEeachFront_.setSize(bDataL_.size());
    setLengthScalesNearTheFront();
}

//note: func to calc the length scale of mesh.
template<class CloudType>
inline Foam::scalar Foam::FrontDynamic<CloudType>::calcLengthScaleOfTheMesh()
{
    /*
    const fvMesh& mesh = this->owner().mesh();
    scalar lengthSum = 0;
  
    scalar counter = 0;
    forAll(mesh.V().field(), cI)
    {
        lengthSum = lengthSum + pow( mesh.V().field()[cI], 1.0/3.0);
        counter = counter + 1.0;
    }
    return lengthSum/mesh.V().field().size();
   */

    const fvMesh& mesh = this->owner().mesh();
    const pointField& ctrs = mesh.cellCentres();
    const scalarField& cv = mesh.V();

    // Integrate variables over cell set
    scalar lengthSum = 0.0;
    scalar counter = 0;
    forAll(ctrs, i)
    {
        label cellI = i;
        scalar volCell = cv[cellI];
        lengthSum += pow( volCell, 1.0/3.0);
        counter += 1.0;
    }

    // Collect across all processors
    reduce(lengthSum, sumOp<scalar>());
    reduce(counter, sumOp<scalar>());

    return lengthSum/counter;
}

//note: func to process the front meshes.
template<class CloudType>
inline void Foam::FrontDynamic<CloudType>::processingThefrontMeshes()
{
    Info << "\n---------> At the time: " << this->owner().db().time().value()
         <<" , Now all front meshes is refining/coarsening  ....... " << endl;

    forAll(bDataL_,bDI)
    {//1
	DynamicList<pointData>& pL = bDataL_[bDI].pL;
        DynamicList<elementInfo>& eL = bDataL_[bDI].eL;

         Info << "      for bubble " << bDI << endl;

         Info << "          Before refining and coarsening, the pL/eL sizes were "
              << pL.size() << ' ' << eL.size() << endl; 

	frontRefining( pL, eL, minEdge_, maxEdge_, maxAspectRatio_);
	frontCoarsening( pL, eL, minEdge_, maxEdge_, maxAspectRatio_);
        adjustingTheMeshDataAfterChanges( pL, eL);

        Info << "          After  refining and coarsening, the pL/eL sizes are "
              << pL.size() << ' ' << eL.size() << endl;//jj
    }//1

    undulationRemovalingOption();
}

//note: func to removal the front undulatuons.
template<class CloudType>
inline void Foam::FrontDynamic<CloudType>::undulationRemovalingOption()
{
    if( URIntervalCounter_ == undulationRemovalingInterval_ && undulationRemovalingOption_ )
    {
        Info << "\n---------> At the time: " << this->owner().db().time().value()
             <<" , Now removaling the front undulations is started  ....... " << endl;

        if( undulationRemovalingType_ == "TSUR3D" )
        {
            forAll(bDataL_,bDI)
            {//1
	        DynamicList<pointData>& pL = bDataL_[bDI].pL;
                DynamicList<elementInfo>& eL = bDataL_[bDI].eL;
	        TSUR3D( 1, pL, eL);
            }//1
        }
        else if( undulationRemovalingType_ == "volumeConservingSmoothingIII" )
        {
            forAll(bDataL_,bDI)
            {//1
	        DynamicList<pointData>& pL = bDataL_[bDI].pL;
                DynamicList<elementInfo>& eL = bDataL_[bDI].eL;
	        volumeConservingSmoothingIII( 1, pL, eL);
            }//1
        }
        else if( undulationRemovalingType_ == "volumeConservingSmoothingIV" || undulationRemovalingType_ == "combined" )
        {
            forAll(bDataL_,bDI)
            {//1
	        DynamicList<pointData>& pL = bDataL_[bDI].pL;
                DynamicList<elementInfo>& eL = bDataL_[bDI].eL;
                //initialBubbleVolume_ = calcBubbleVolume(eL);
	        volumeConservingSmoothingIV( 1, pL, eL);
            }//1
        }
        else
        {
            FatalErrorIn ("Foam::FrontDynamic<CloudType>::undulationRemovalingOption()")
            << "you should choose undulationRemovalingType as TSUR3D or volumeConservingSmoothingIII or volumeConservingSmoothingIV" << nl
            << abort(FatalError);
        }
        URIntervalCounter_ = 0;
    }
    URIntervalCounter_++;         
}

//note: func to removal the front undulatuons.
template<class CloudType>
inline void Foam::FrontDynamic<CloudType>::adjustingTheMeshDataAfterChanges(DynamicList<pointData>& pL, DynamicList<elementInfo>& eL)
{
    DynamicList<pointData> pLTemp = pL;
    DynamicList<elementInfo> eLTemp = eL;

    //DynamicList<pointData> pLTemp ;
    //DynamicList<elementInfo> eLTemp;
    //pLTemp.setCapacity(pL.size());
    //eLTemp.setCapacity(eL.size());
    //pLTemp = pL;
    //eLTemp = eL;

    //locInPointDataList_.resize( pL.size());
    //locInElementInfoList_.resize( eL.size());
    locInPointDataList_.setSize( pL.size());
    locInElementInfoList_.setSize( eL.size());

    //resizing and renumbering the point and element lists.
    //pL.clearStorage();
    //eL.clearStorage();
    pL.clear();//JJ
    eL.clear();//JJ

    label ptCount = 0;
    forAll(pLTemp, I)
    {//2
        pointData& ptI = pLTemp[I];
        locInPointDataList_[ptI.currentIndex] = -1;
        if( ptI.keepPoint == true)
        {
            locInPointDataList_[ptI.currentIndex] = ptCount;
            ptI.currentIndex = ptCount;
            pL.append(ptI);
            ptCount = ptCount + 1;
        }
    }//2

    label elCount = 0;
    forAll(eLTemp, II)
    {//3
        elementInfo& elI = eLTemp[II];
        locInElementInfoList_[elI.currentIndex] = -1;
        if( elI.keepElement == true)
        {
            locInElementInfoList_[elI.currentIndex] = elCount; 
            elI.currentIndex = elCount;
            eL.append(elI);
            elCount = elCount + 1;
        }
    }//3

    //Info << "          After  refining and coarsening, the pL/eL sizes are " << ptCount << ' ' << elCount << "\n" << endl;//JJ

    pLTemp.clearStorage();
    eLTemp.clearStorage();

    forAll(eL, II)
    {//4
        elementInfo& elI = eL[II];

        label& p0 = elI.pointIndex[0];
        label& p1 = elI.pointIndex[1];
        label& p2 = elI.pointIndex[2];

        elI.pointIndex[0] = locInPointDataList_[p0];
        elI.pointIndex[1] = locInPointDataList_[p1];
        elI.pointIndex[2] = locInPointDataList_[p2];

        elI.points[0] = pL[p0].currentPoint;
        elI.points[1] = pL[p1].currentPoint;
        elI.points[2] = pL[p2].currentPoint;

        label e0 = elI.elementIndex[0];
        label e1 = elI.elementIndex[1];
        label e2 = elI.elementIndex[2];

        elI.elementIndex[0] = locInElementInfoList_[e0];
        elI.elementIndex[1] = locInElementInfoList_[e1];
        elI.elementIndex[2] = locInElementInfoList_[e2];

        pL[p0].elOwner = elI.currentIndex;
        pL[p1].elOwner = elI.currentIndex;
        pL[p2].elOwner = elI.currentIndex;
    }//4

    ptConnectedPoints(pL,eL);

    locInPointDataList_.clear();
    locInElementInfoList_.clear();           
}


//note: func for adjusting the important vars for indicators.
template<class CloudType>
inline void Foam::FrontDynamic<CloudType>::adjustingTheImportantVarsForIndicators()
{
    p_[0] = 0;
    p_[1] = 1;
    p_[2] = 0;
    p_[3] = 0;
    p_[4] = 1;
    p_[5] = 2;
    p_[6] = 1;
    p_[7] = 1;
    p_[8] = 0;
    p_[9] = 1;
    p_[10] = 0;
    p_[11] = 0;
    p_[12] = 0;
    p_[13] = 1;
    p_[14] = 0;
    p_[15] = 0;
 
    q_[0] = 0;
    q_[1] = 0;
    q_[2] = 1;
    q_[3] = 0;
    q_[4] = 0;
    q_[5] = 0;
    q_[6] = 1;
    q_[7] = 0;
    q_[8] = 1;
    q_[9] = 1;
    q_[10] = 2;
    q_[11] = 1;
    q_[12] = 0;
    q_[13] = 0;
    q_[14] = 1;
    q_[15] = 0;

    r_[0] = 0;
    r_[1] = 0;
    r_[2] = 0;
    r_[3] = 1;
    r_[4] = 0;
    r_[5] = 0;
    r_[6] = 0;
    r_[7] = 1;
    r_[8] = 0;
    r_[9] = 0;
    r_[10] = 0;
    r_[11] = 1;
    r_[12] = 1;
    r_[13] = 1;
    r_[14] = 1;
    r_[15] = 2;

    //mComponent_.source()[0] = 1.0;
    //mComponent_.source()[1] = 0.0;
    //mComponent_.source()[2] = 0.0;
    //mComponent_.source()[3] = 0.0;
}

//note: func for constructing the indicators based on CPT.
template<class CloudType>
inline void Foam::FrontDynamic<CloudType>::constructingTheIndicatorsBasedOnCPT()
{
    Info << "\n---------> At the time: " << this->owner().db().time().value()
         <<" , Now densityIndicator_ and viscosityIndicator are evaluated based on CPT method  .......\n " << endl;

    const fvMesh& mesh = this->owner().mesh();
    const pointField& ctrs = mesh.cellCentres();

  //searchAllCellsNearEeachFront();

    forAll(bDataL_, bDI)
    {//10
        DynamicList<elementInfo>& eL = bDataL_[bDI].eL;
        DynamicList<label>& ctrsLabel = allCellsInEeachMasket_[bDI];
      //DynamicList<label>& ctrsLabel = allCellsNearEeachFront_[bDI];
        h_ = lengthScaleNearTheFront_[bDI];

        forAll(ctrsLabel,I)
        {//5
            const label& cI = ctrsLabel[I];
 
            CPTData& CPT = CPTList_[cI];

            vector minDistance; 
            label targetElement;
            minDistance.x() = GREAT;   
            minDistance.y() = GREAT;   
            minDistance.z() = GREAT;   

	    forAll(eL, elementI)
	    {//4
	        const elementInfo& elInfo = eL[elementI];
		const point& pos = elInfo.centrePosInDomain;
		vector distance = ctrs[cI] - pos;
                if( mag(distance) <= mag(minDistance) )
                {
                    targetElement = elementI;
                    minDistance = distance;
                }
            }//4
            
            // calculating the accurate closest point in each element
            //calculatingTheExactClosestPointForTheCurrentElement(eL[targetElement], minDistance);

            //calculating the density and viscosity induced by each bubble.
            CPT.inducedDensity[bDI]= LiuHeviside(bDataL_[bDI].outerFluidDensity, bDataL_[bDI].density,
                                                  minDistance, eL[targetElement].elementSurfaceArea, 2.*h_);
            CPT.inducedViscosity[bDI]= LiuHeviside(bDataL_[bDI].outerFluidViscosity, bDataL_[bDI].viscosity,
                                                  minDistance, eL[targetElement].elementSurfaceArea, 2.*h_);
        }//5
    }//10

    // constructing final indicators
    forAll(ctrs, cI)
    {//5
	CPTData& CPT = CPTList_[cI];
		
	scalar finalDensity = 0.0;
	scalar finalViscosity = 0.0;
	scalar densityDiffMax = 0.0;
	scalar viscosityDiffMax = 0.0;

	for(int bDI = 0; bDI < bDataL_.size() ; bDI++)
	{//0
	    scalar densityDiff = mag(CPT.inducedDensity[bDI] -bDataL_[bDI].outerFluidDensity);
	    scalar viscosityDiff = mag(CPT.inducedViscosity[bDI] -bDataL_[bDI].outerFluidViscosity); 

	    if( densityDiff >= densityDiffMax )
	    {
		densityDiffMax = densityDiff;
	        finalDensity = CPT.inducedDensity[bDI];
	    }
	    if( viscosityDiff >= viscosityDiffMax )
	    {
		viscosityDiffMax = viscosityDiff;
	        finalViscosity = CPT.inducedViscosity[bDI];
	    }
	}//0

	//CPT.finalDensity = finalDensity;
	//CPT.finalViscosity = finalViscosity;
	densityIndicator_.internalField()[cI] = finalDensity;
	viscosityIndicator_.internalField()[cI] = finalViscosity;
    }//5
}

//note: func for calculating the exact closest point at the current element.
template<class CloudType>
inline void Foam::FrontDynamic<CloudType>::calculatingTheExactClosestPointForTheCurrentElement(elementInfo elInfo, vector& minDistance)
{
    point p1 = elInfo.points[0];
    point p2 = elInfo.points[1];
    point p3 = elInfo.points[2];

    point mainElementCentre = ( p1 + p2 + p3)/3.0;

    point gridPos = mainElementCentre + minDistance;

    for(int i = 1; i <= numberOfSearchToObtainTheExactCPForEachElement_ ; i++) 
    {
        mainElementCentre = ( p1 + p2 + p3)/3.0;
        point newElementICentre = ( p1 + p2 + mainElementCentre)/3.0;
        point newElementIICentre = ( p2 + p3 + mainElementCentre)/3.0;
        point newElementIIICentre = ( p3 + p1 + mainElementCentre)/3.0;

        scalar d1 = mag( gridPos - newElementICentre);
        scalar d2 = mag( gridPos - newElementIICentre); 
        scalar d3 = mag( gridPos - newElementIIICentre);

        if( d1 <= mag(minDistance) && d1 <= d2 && d1 <= d3 )
        {
            p1 = p1;
            p2 = p2;
            p3 = mainElementCentre;
            minDistance = gridPos - newElementICentre; 
        }
        else if( d2 <= mag(minDistance) && d2 <= d1 && d2 <= d3 )
        {
            p1 = p2;
            p2 = p3;
            p3 = mainElementCentre;
            minDistance = gridPos - newElementIICentre; 
        }
        else if( d3 <= mag(minDistance) && d3 <= d2 && d3 <= d1 )
        {
            p1 = p3;
            p2 = p1;
            p3 = mainElementCentre;
            minDistance = gridPos - newElementIIICentre; 
        }
        else
        {
            break;
        }
    }
}

//note: func for calculating the exact closest point at the current element.
template<class CloudType>
inline Foam::scalar Foam::FrontDynamic<CloudType>::LiuHeviside( scalar outerFluidProperty, scalar property, vector minDistance, vector elementSurfaceArea, scalar gama)
{
    // indicators are built with the aid of a Heaviside function (Liu et al., 2000)

    if ( mag(minDistance) <= gama && (minDistance & elementSurfaceArea) > 0.0 )
    {
        scalar phi = mag(minDistance);
        scalar inducedProperty =  0.5 * (outerFluidProperty - property)
                                      * (1.0 + phi/gama + sin(mathematical::pi * phi/gama)/mathematical::pi) + property;
        return inducedProperty;
    }
    if ( mag(minDistance) <= gama && (minDistance & elementSurfaceArea) < 0.0 )
    {
        scalar phi = - mag(minDistance);
        scalar inducedProperty =  0.5 * (outerFluidProperty - property)
                                         * (1.0 + phi/gama + sin(mathematical::pi * phi/gama)/mathematical::pi) + property;
        return inducedProperty;
    } 
    if ( mag(minDistance) > gama && (minDistance & elementSurfaceArea) > 0.0 )
    {
        return outerFluidProperty;
    }
    if ( mag(minDistance) > gama && (minDistance & elementSurfaceArea) < 0.0 )
    {
        return property;
    }
}

//note: func for search all cells existed in eeach masket.
template<class CloudType>
inline void Foam::FrontDynamic<CloudType>::searchAllCellsInEeachMasket()
{
    const fvMesh& mesh = this->owner().mesh();
    const pointField& ctrs = mesh.cellCentres();
    const scalarField& cv = mesh.V();

    initialAllocatingTheCPTList();

    forAll(bDataL_,bDI)
    {//10
	const DynamicList<pointData>& pL = bDataL_[bDI].pL;
        DynamicList<label>& L = allCellsInEeachMasket_[bDI];
        L.clearStorage();

        h_ = lengthScaleNearTheFront_[bDI];
        scalar& lengthScale = lengthScaleNearTheFront_[bDI];

        scalar xMin = +100000000.;
        scalar xMax = -100000000.;
        scalar yMin = +100000000.;
        scalar yMax = -100000000.;
        scalar zMin = +100000000.;
        scalar zMax = -100000000.;

        for( int i=0; i < pL.size(); i++)
        {
            // search for xMax and xMin
            if( pL[i].currentPoint.x() <= xMin )
            {
                xMin = pL[i].currentPoint.x();
            }
            if( pL[i].currentPoint.x() >= xMax )
            {
                xMax = pL[i].currentPoint.x();
            }

            // search for yMax and yMin
            if( pL[i].currentPoint.y() <= yMin )
            {
                yMin = pL[i].currentPoint.y();
            }
            if( pL[i].currentPoint.y() >= yMax )
            {
                yMax = pL[i].currentPoint.y();
            }

            // search for zMax and zMin
            if( pL[i].currentPoint.z() <= zMin )
            {
                zMin = pL[i].currentPoint.z();
            }
            if( pL[i].currentPoint.z() >= zMax )
            {
                zMax = pL[i].currentPoint.z();
            }
        }

        Info << " xMin/xMax " << xMin << ' ' << xMax << endl;
        Info << " yMin/yMax " << yMin << ' ' << yMax << endl;
        Info << " zMin/zMax " << zMin << ' ' << zMax << endl;

        xMin = xMin - 3.0 * h_;
        xMax = xMax + 3.0 * h_;
        yMin = yMin - 3.0 * h_;
        yMax = yMax + 3.0 * h_;
        zMin = zMin - 3.0 * h_;
        zMax = zMax + 3.0 * h_;

        // calc the length scale near the front.
        scalar lengthSum = 0.0;
        scalar counter = 0;

        forAll(ctrs, cI)
        {//2
            CPTData& CPT = CPTList_[cI];
            if( (xMin <= ctrs[cI].x() && ctrs[cI].x() <= xMax) &&
                (yMin <= ctrs[cI].y() && ctrs[cI].y() <= yMax) &&
                (zMin <= ctrs[cI].z() && ctrs[cI].z() <= zMax)  )
            {
                //Indicator_.internalField()[cI] = 0.0;
                CPT.inducedDensity[bDI] = bDataL_[bDI].density;
                CPT.inducedViscosity[bDI] = bDataL_[bDI].viscosity;
                L.append(cI);

                // calc the length scale near the front.
                scalar volCell = cv[cI];
                lengthSum += pow( volCell, 1.0/3.0);
                counter += 1.0; 
            }
            //Indicator_.internalField()[cI] = 1.0;
            CPT.inducedDensity[bDI] = bDataL_[bDI].outerFluidDensity;
            CPT.inducedViscosity[bDI] = bDataL_[bDI].outerFluidViscosity;
        }//2

        // calc the length scale near the front.
        reduce(lengthSum, sumOp<scalar>());
        reduce(counter, sumOp<scalar>());
        lengthScale = lengthSum/counter;
    }//10
}

//note: func for searchAllCellsNearEeachFront.
template<class CloudType>
inline void Foam::FrontDynamic<CloudType>::searchAllCellsNearEeachFront()
{
/*
    const fvMesh& mesh = this->owner().mesh();
    const pointField& ctrs = mesh.cellCentres();

    forAll(bDataL_,bDI)
    {//10
        DynamicList<label>& L = allCellsInEeachMasket_[bDI];
        DynamicList<label>& LL = allCellsNearEeachFront_[bDI];
        LL.clearStorage();

        forAll(L, LI)
        {//15
            const label& cI = L[LI];
            CPTData& CPT = CPTList_[cI];
            if( cellsToRefine_.internalField()[cI] == -1. )
            {
                LL.append(cI);
            }
            else if( CPT.inducedDensity[bDI] != bDataL_[bDI].density && cellsToRefine_.internalField()[cI] == 0.0 )
            {
                CPT.inducedDensity[bDI] = bDataL_[bDI].outerFluidViscosity;
                CPT.inducedViscosity[bDI] = bDataL_[bDI].outerFluidViscosity;
            }
        }//15
    }//10
*/
}

//note: func for communications of the front and Eulerian grid.
template<class CloudType>
inline void Foam::FrontDynamic<CloudType>::communicationsOfTheFrontAndEulerianGrid()
{
    Info << "\n---------> At the time: " << this->owner().db().time().value()
         <<" , Now all communications of the front and Eulerian grid are doing  .......\n " << endl;

    const fvMesh& mesh = this->owner().mesh();
    const pointField& ctrs = mesh.cellCentres();

    this->owner().sTensionForceFromFT().field() = vector::zero;
    this->owner().pressureJumpAtTheInterfaceFromFT().field() = vector::zero;
    cellsToRefine_.field() = 0.0;
    searchAllCellsInEeachMasket();

    simpleMatrix<scalar> mComponent(4);
    mComponent.source()[0] = 1.0;
    mComponent.source()[1] = 0.0;
    mComponent.source()[2] = 0.0;
    mComponent.source()[3] = 0.0;

    volSymmTensorField shearStress = 2*viscosityIndicator_
                                      *symm( fvc::grad(this->owner().U()) );

    volVectorField densityIGrad
    (
        IOobject
        (
            "densityIGrad",
            this->owner().db().time().timeName(),
            mesh
        ),
        mesh,
        dimensionedVector
        (
            "densityIGrad",
            densityIndicator_.dimensions()/dimLength,
            vector::zero
        ),
        zeroGradientFvPatchVectorField::typeName
    );
    densityIGrad.internalField() = vector::zero;

    volVectorField viscosityIGrad
    (
        IOobject
        (
            "viscosityIGrad",
            this->owner().db().time().timeName(),
            mesh
        ),
        mesh,
        dimensionedVector
        (
            "viscosityIGrad",
            viscosityIndicator_.dimensions()/dimLength,
            vector::zero
        ),
        zeroGradientFvPatchVectorField::typeName
    );
    viscosityIGrad.internalField() = vector::zero;

    volVectorField IndicatorGrad
    (
        IOobject
        (
            "IndicatorGrad",
            this->owner().db().time().timeName(),
            mesh
        ),
        mesh,
        dimensionedVector
        (
            "IndicatorGrad",
            Indicator_.dimensions()/dimLength,
            vector::zero
        ),
        zeroGradientFvPatchVectorField::typeName
    );
    IndicatorGrad.internalField() = vector::zero;

    DynamicList<label> cellsInSphereDomain;

    if (!Pstream::parRun())
    {//5
        forAll(bDataL_,bDI)
        {//10
	    DynamicList<elementInfo>& eL = bDataL_[bDI].eL;
            DynamicList<label>& ctrsLabel = allCellsInEeachMasket_[bDI];
            h_ = lengthScaleNearTheFront_[bDI];

	    forAll(eL, elementI)
	    {//4
		elementInfo& elInfo = eL[elementI];   
		const point& pos = elInfo.centrePosInDomain;

                forAll(ctrsLabel, I)
                {
                    const label& cI = ctrsLabel[I];
                    scalar distance = mag (ctrs[cI] - pos);
                    if ( distance <= 2.0 * h_ )
                    {
                        cellsInSphereDomain.append(cI); 
                    }
                }
                
		for(int i=0; i < 4 ; i++)
		{//1
		    for(int j=0; j < 4 ; j++)
		    {//2
			 label pqr = 4*i + j;
			 mComponent[i][j] = 0.0;
			 forAll(cellsInSphereDomain, cI)
			 {
			     label cellI = cellsInSphereDomain[cI];
			     vector del = ( ctrs[cellI] - pos)/(h_ + ROOTVSMALL);
                             scalar m = mPQR(del, h_, p_[pqr], q_[pqr], r_[pqr]);
			     mComponent[i][j] = mComponent[i][j] + m;
			 }

		    }//2
		}//1

		beta_ = mComponent.solve();
                elInfo.pressure = 0.0;
                elInfo.shearStress = symmTensor::zero;
		forAll(cellsInSphereDomain, cI)
		{
		    label cellI = cellsInSphereDomain[cI];
		    scalar dVolume = mesh.V()[cellI];
		    vector del = ( ctrs[cellI] - pos)/(h_ + ROOTVSMALL);
                    scalar wBar = windowFuncBar(del, h_, beta_);
		    densityIGrad.internalField()[cellI] = densityIGrad.internalField()[cellI] 
                                                   + (bDataL_[bDI].outerFluidDensity-bDataL_[bDI].density)  
                                                   * wBar * elInfo.elementSurfaceArea/dVolume;
		    viscosityIGrad.internalField()[cellI] = viscosityIGrad.internalField()[cellI]
                                                   + (bDataL_[bDI].outerFluidViscosity-bDataL_[bDI].viscosity)
                                                   * wBar * elInfo.elementSurfaceArea/dVolume;
		    IndicatorGrad.internalField()[cellI] = IndicatorGrad.internalField()[cellI]
                                                   + (1.0-0.0)
                                                   * wBar * elInfo.elementSurfaceArea/dVolume;

                    elInfo.pressure = elInfo.pressure + wBar * (this->owner().p().internalField()[cellI]);
                    elInfo.shearStress = elInfo.shearStress + wBar * (shearStress.internalField()[cellI]);

                    this->owner().sTensionForceFromFT()[cellI] = this->owner().sTensionForceFromFT()[cellI]
                                     + wBar * elInfo.averageSurfaceTension * mag(elInfo.elementSurfaceArea)/dVolume;

	            if( pressureJumpAtTheInterfaceOption_ == true)
	            {
                        this->owner().pressureJumpAtTheInterfaceFromFT()[cellI] = this->owner().pressureJumpAtTheInterfaceFromFT()[cellI]
                                     + wBar * bDataL_[bDI].pressureJumpAtTheInterface * mag(elInfo.elementSurfaceArea)/dVolume;
	            }

		    cellsToRefine_.internalField()[cellI] = -1.0;
		}
		cellsInSphereDomain.clear();
	    }//4
        }//10
    }//5

    if (Pstream::parRun())
    {//5
        List<scalar> mComponentL(16);

        forAll(bDataL_,bDI)
        {//10
	    DynamicList<elementInfo>& eL = bDataL_[bDI].eL;
            DynamicList<label>& ctrsLabel = allCellsInEeachMasket_[bDI];
            h_ = lengthScaleNearTheFront_[bDI];       

	    forAll(eL, elementI)
	    {//4
		elementInfo& elInfo = eL[elementI];   
		const point& pos = elInfo.centrePosInDomain;

                forAll(ctrsLabel, I)
                {
                    const label& cI = ctrsLabel[I];
                    scalar distance = mag (ctrs[cI] - pos);
                    if ( distance <= 2.0 * h_ )
                    {
                        cellsInSphereDomain.append(cI); 
                    }
                }

		for(int i=0; i < 4 ; i++)
		{//1
		    for(int j=0; j < 4 ; j++)
		    {//2
			 label pqr = 4*i + j;
			 mComponentL[pqr] = 0.0;
			 mComponent[i][j] = 0.0;
			 forAll(cellsInSphereDomain, cI)
			 {
			     label cellI = cellsInSphereDomain[cI];
			     vector del = ( ctrs[cellI] - pos)/(h_ + ROOTVSMALL);
                             scalar m = mPQR(del, h_, p_[pqr], q_[pqr], r_[pqr]);
			     mComponentL[pqr] = mComponentL[pqr] + m;
			 }
		    }//2
		}//1

		// gather/scatter for mComponentL
		List<List<scalar> > mComponentLAll(Pstream::nProcs());
		mComponentLAll[Pstream::myProcNo()] = mComponentL;
		Pstream::gatherList(mComponentLAll);
		Pstream::scatterList(mComponentLAll);
		forAll(mComponentLAll,mCLAI)
	        {
	            List<scalar> tempMComponentL = mComponentLAll[mCLAI];
	            for(int i=0; i < 4 ; i++)
	            {//1
	                for(int j=0; j < 4 ; j++)
		        {//2
	                    label pqr = 4*i + j;
		            mComponent[i][j] = mComponent[i][j] + tempMComponentL[pqr];
		        }//2
	            }//1
	        }

		beta_ = mComponent.solve();
                elInfo.pressure = 0.0;
                elInfo.shearStress = symmTensor::zero;
		forAll(cellsInSphereDomain, cI)
		{
		    label cellI = cellsInSphereDomain[cI];
		    scalar dVolume = mesh.V()[cellI];
		    vector del = ( ctrs[cellI] - pos)/(h_ + ROOTVSMALL);
                    scalar wBar = windowFuncBar(del, h_, beta_);
		    densityIGrad.internalField()[cellI] = densityIGrad.internalField()[cellI] 
                                                   + (bDataL_[bDI].outerFluidDensity-bDataL_[bDI].density)  
                                                   * wBar * elInfo.elementSurfaceArea/dVolume;
		    viscosityIGrad.internalField()[cellI] = viscosityIGrad.internalField()[cellI]
                                                   + (bDataL_[bDI].outerFluidViscosity-bDataL_[bDI].viscosity)
                                                   * wBar * elInfo.elementSurfaceArea/dVolume;
		    IndicatorGrad.internalField()[cellI] = IndicatorGrad.internalField()[cellI]
                                                   + (1.0-0.0)
                                                   * wBar * elInfo.elementSurfaceArea/dVolume;

                    elInfo.pressure = elInfo.pressure + wBar * (this->owner().p().internalField()[cellI]);
                    elInfo.shearStress = elInfo.shearStress + wBar * (shearStress.internalField()[cellI]);

                    this->owner().sTensionForceFromFT()[cellI] = this->owner().sTensionForceFromFT()[cellI]
                                     + wBar * elInfo.averageSurfaceTension * mag(elInfo.elementSurfaceArea)/dVolume;

	            if( pressureJumpAtTheInterfaceOption_ == true)
	            {
                        this->owner().pressureJumpAtTheInterfaceFromFT()[cellI] = this->owner().pressureJumpAtTheInterfaceFromFT()[cellI]
                                     + wBar * bDataL_[bDI].pressureJumpAtTheInterface * mag(elInfo.elementSurfaceArea)/dVolume;
	            }

		    cellsToRefine_.internalField()[cellI] = -1.0;
		}
		cellsInSphereDomain.clear();
	    }//4
        }//10
    }//5

    if ( constructingTheIndicatorsOption_ == "basedOnPoissonEquation" )
    { 
        if ( twoFluidFlow_ == true )
        {
            IndicatorGrad.correctBoundaryConditions();

            fvScalarMatrix IndicatorEqn
            (
                fvm::laplacian(Indicator_)
	      ==
                fvc::div(IndicatorGrad)
            );
            IndicatorEqn.solve();

            forAll(ctrs, cI)
            {
	        if( Indicator_.internalField()[cI] > 0.998 )
	        {
	            Indicator_.internalField()[cI] = 1.0;
	        }
            }
            Indicator_ = max(Indicator_, 0.0);
            Indicator_ = min(Indicator_, 1.0);
	    densityIndicator_ = (1.0 - Indicator_) * bDataL_[0].density + Indicator_ * bDataL_[0].outerFluidDensity;
	    viscosityIndicator_ = (1.0 - Indicator_) * bDataL_[0].viscosity + Indicator_ * bDataL_[0].outerFluidViscosity;
        }
        else if ( twoFluidFlow_ == false )
        {
            densityIGrad.correctBoundaryConditions();
            viscosityIGrad.correctBoundaryConditions();

            fvScalarMatrix densityIEqn
            (
                fvm::laplacian(densityIndicator_)
	      ==
                fvc::div(densityIGrad)
            );
            densityIEqn.solve();

            fvScalarMatrix viscosityIEqn
            (
                fvm::laplacian(viscosityIndicator_)
	      ==
                fvc::div(viscosityIGrad)
            );
            viscosityIEqn.solve();
        }
    }
    else if ( constructingTheIndicatorsOption_ == "basedOnCPT" )
    {
        constructingTheIndicatorsBasedOnCPT();
    }

    for(int i = 1; i <= numberOfIndicatorFiltering_ ; i++)
    {
        densityIndicator_ = filteringWeight_ * (fvc::average(fvc::interpolate(densityIndicator_))) 
	                                     + (1.0 - filteringWeight_) * densityIndicator_;
	viscosityIndicator_ = filteringWeight_ * (fvc::average(fvc::interpolate(viscosityIndicator_))) 
	                                     + (1.0 - filteringWeight_) * viscosityIndicator_;
    } 
    this->owner().densityFromFT() = densityIndicator_;
    this->owner().viscosityFromFT() = viscosityIndicator_;
}

//note: func for surface Tension Force Distribution.
template<class CloudType>
inline void Foam::FrontDynamic<CloudType>::surfaceTensionForceDistribution()
{
    forAll(bDataL_,bDI)
    {//10
        vector sumSurfaceTension = vector::zero;
	DynamicList<elementInfo>& eL = bDataL_[bDI].eL;

	forAll(eL, elementI)
	{//4
	    elementInfo& elInfo = eL[elementI];
            elInfo.elementSurfaceArea = calcElementSurfaceArea(elInfo);   
        }//4

	forAll(eL, elementI)
	{//5
	    elementInfo& elInfo = eL[elementI];   

	    if( surfaceTensionForceDistributionOption_ == "basedOnElementNeighbours")
	    {
	        label p1I = elInfo.pointIndex[0];
	        label p2I = elInfo.pointIndex[1];
	        label p3I = elInfo.pointIndex[2];

	        point p1 = bDataL_[bDI].pL[p1I].currentPoint;
	        point p2 = bDataL_[bDI].pL[p2I].currentPoint;
	        point p3 = bDataL_[bDI].pL[p3I].currentPoint;

	        point t1 = p2 - p1;
	        point t2 = p3 - p2;
	        point t3 = p1 - p3;

	        label e1I = elInfo.elementIndex[0];
	        label e2I = elInfo.elementIndex[1];
	        label e3I = elInfo.elementIndex[2];

	        point ne1 = eL[e1I].elementSurfaceArea/mag(eL[e1I].elementSurfaceArea);
	        point ne2 = eL[e2I].elementSurfaceArea/mag(eL[e2I].elementSurfaceArea);
	        point ne3 = eL[e3I].elementSurfaceArea/mag(eL[e3I].elementSurfaceArea);
                
                elInfo.averageSurfaceTension = 0.5 * bDataL_[bDI].surfaceTensionCoeff * ( (t1 ^ ne1) + (t2 ^ ne2) + (t3 ^ ne3) );
	    }
	    else if( surfaceTensionForceDistributionOption_ == "basedOnPointNeighbours")
	    {
	        label p1I = elInfo.pointIndex[0];
	        label p2I = elInfo.pointIndex[1];
	        label p3I = elInfo.pointIndex[2];
                elInfo.averageSurfaceTension = (  calcSurfaceTensionForceAtpoint( bDataL_[bDI].pL[p1I])
                                                + calcSurfaceTensionForceAtpoint( bDataL_[bDI].pL[p2I])
                                                + calcSurfaceTensionForceAtpoint( bDataL_[bDI].pL[p3I]) )/3.0;
                elInfo.averageSurfaceTension = elInfo.averageSurfaceTension * bDataL_[bDI].surfaceTensionCoeff;
	    }

            sumSurfaceTension = sumSurfaceTension + elInfo.averageSurfaceTension;
	}//5

        bDataL_[bDI].pressureJumpAtTheInterface = sumSurfaceTension/calcBubbleSurface(eL);
    }//10
}

//-note: func for catching the required fields
template<class CloudType>
inline void Foam::FrontDynamic<CloudType>::cacheFields()
{
    //tmp<volScalarField> tp_rgh = p_rgh();
    //p_rghPtr_ = tp_rgh.ptr();

    tmp<volScalarField> tp = p();
    pPtr_ = tp.ptr();

    tmp<volVectorField> tU = U();
    UPtr_ = tU.ptr();
}

//-note: func for accessing the p_rgh field.
template<class CloudType>
inline Foam::tmp<Foam::volScalarField> Foam::FrontDynamic<CloudType>::p_rgh() const
{
    const fvMesh& mesh = this->owner().mesh();

    const objectRegistry& obr = mesh;

    if (obr.foundObject<volScalarField>("p_rgh"))
    {
		const volScalarField& p_rghc =
            obr.lookupObject<volScalarField>("p_rgh");

        return p_rghc;
    }
    else
    {
        FatalErrorIn
        (
            "inline Foam::tmp<Foam::volScalarField> Foam::FrontDynamic<CloudType>::p_rgh() const"
        )
            << "carrier phase fields not found in mesh database" << nl
            << "Database objects include: " << obr.sortedToc()
            << abort(FatalError);

        return tmp<volScalarField>(NULL);
    }
}

//-note: func for accessing the p = p_rgh+rgh field.
template<class CloudType>
inline Foam::tmp<Foam::volScalarField> Foam::FrontDynamic<CloudType>::p() const
{
    const fvMesh& mesh = this->owner().mesh();

    const objectRegistry& obr = mesh;

    if (obr.foundObject<volScalarField>("p"))
    {
		const volScalarField& pc =
            obr.lookupObject<volScalarField>("p");

        return pc;
    }
    else
    {
        FatalErrorIn
        (
            "inline Foam::tmp<Foam::volScalarField> Foam::FrontDynamic<CloudType>::p() const"
        )
            << "carrier phase fields not found in mesh database" << nl
            << "Database objects include: " << obr.sortedToc()
            << abort(FatalError);

        return tmp<volScalarField>(NULL);
    }
}

//-note: func for accessing the U field.
template<class CloudType>
inline Foam::tmp<Foam::volVectorField> Foam::FrontDynamic<CloudType>::U() const
{
    const fvMesh& mesh = this->owner().mesh();

    const objectRegistry& obr = mesh;

    if (obr.foundObject<volVectorField>("U"))
    {
		const volVectorField& Uc =
            obr.lookupObject<volVectorField>("U");

        return Uc;
    }
    else
    {
        FatalErrorIn
        (
            "inline Foam::tmp<Foam::volVectorField> Foam::FrontDynamic<CloudType>::U() const"
        )
            << "carrier phase fields not found in mesh database" << nl
            << "Database objects include: " << obr.sortedToc()
            << abort(FatalError);

        return tmp<volVectorField>(NULL);
    }
}

//note: func for bubble post-processing.
template<class CloudType>
inline void Foam::FrontDynamic<CloudType>::bubblePostProcessing()
{
    
    Info << "\n----------> At the time: " << this->owner().db().time().value()
         <<" , Now post-processing for all bubbles is doing  ....... " << endl;

    scalar time = this->owner().db().time().value();
    scalar deltaTime = this->owner().db().time().deltaTValue();

    forAll(bDataL_,bDI)
    {//1
	DynamicList<pointData> pL = bDataL_[bDI].pL;
        DynamicList<elementInfo> eL = bDataL_[bDI].eL;

   	bDataL_[bDI].volume = calcBubbleVolume(eL);
   	bDataL_[bDI].surfaceArea = calcBubbleSurface(eL);
        //bDataL_[bDI].sphereBubbleDiameter = pow(6.0*bDataL_[bDI].volume/mathematical::pi, 1./3.);
   	bDataL_[bDI].sphericity = mathematical::pi * sqr( bDataL_[bDI].sphereBubbleDiameter)/bDataL_[bDI].surfaceArea;
   	bDataL_[bDI].mass = bDataL_[bDI].density * bDataL_[bDI].volume;
   	bDataL_[bDI].gravityForce = this->owner().g().value() * bDataL_[bDI].mass;
   	bDataL_[bDI].buoyantForce = this->owner().g().value() * bDataL_[bDI].outerFluidDensity * bDataL_[bDI].volume;

   	bDataL_[bDI].pressureForce = vector::zero;
   	bDataL_[bDI].shearStressForce = vector::zero;
   	bDataL_[bDI].centre = vector::zero;
        scalar sumSurfaceAreaX = 0.0;
        scalar sumSurfaceAreaY = 0.0;
        scalar sumSurfaceAreaZ = 0.0;
   	forAll(eL, elementI)
	{//4
	    elementInfo& elInfo = eL[elementI];

            vector elementSurfaceArea = calcElementSurfaceArea(elInfo);
            sumSurfaceAreaX = sumSurfaceAreaX + mag(elementSurfaceArea.x());
            sumSurfaceAreaY = sumSurfaceAreaY + mag(elementSurfaceArea.y());
            sumSurfaceAreaZ = sumSurfaceAreaZ + mag(elementSurfaceArea.z());

            bDataL_[bDI].pressureForce = bDataL_[bDI].pressureForce - elInfo.pressure * elInfo.elementSurfaceArea;

            bDataL_[bDI].shearStressForce = bDataL_[bDI].shearStressForce + (elInfo.shearStress & elInfo.elementSurfaceArea);

            vector elementCentre = ( pL[elInfo.pointIndex[0]].currentPoint
                                             + pL[elInfo.pointIndex[1]].currentPoint
                                                        + pL[elInfo.pointIndex[2]].currentPoint )/3.0;

            bDataL_[bDI].centre.x() = bDataL_[bDI].centre.x() + elementCentre.x() * mag(elementSurfaceArea.x());
            bDataL_[bDI].centre.y() = bDataL_[bDI].centre.y() + elementCentre.y() * mag(elementSurfaceArea.y());
            bDataL_[bDI].centre.z() = bDataL_[bDI].centre.z() + elementCentre.z() * mag(elementSurfaceArea.z());
        }
        bDataL_[bDI].centre.x() = bDataL_[bDI].centre.x()/sumSurfaceAreaX;
        bDataL_[bDI].centre.y() = bDataL_[bDI].centre.y()/sumSurfaceAreaY;
        bDataL_[bDI].centre.z() = bDataL_[bDI].centre.z()/sumSurfaceAreaZ;

        // calculating the flow direction:
	if( periodicalOption_ == "rotationalPeriodicalInCurvedDuct" || periodicalOption_ == "combinedRotationalPeriodicalInCurvedDuct" ||
            periodicalOption_ == "combinedTranslationalPeriodicalInCurvedDuct" || periodicalOption_ == "translationalPeriodicalInCurvedDuct" )
	{
	    rCentre_.z() =  bDataL_[bDI].centre.z();
	    vector rStar = bDataL_[bDI].centre - rCentre_;

            normalToflowDI_ = rStar/mag(rStar);

	    flowD_.x() = -rStar.y();
	    flowD_.y() =  rStar.x();
	    flowD_.z() =  0;
	    flowD_ = flowD_/mag(flowD_);

	    normalToflowDII_.x() = 0;
	    normalToflowDII_.y() = 0;
	    normalToflowDII_.z() = 1;
	} 
        else if( periodicalOption_ == "straightDuctPeriodical" || periodicalOption_ == "none")
        {
            // for none option, flowD is parallel to z axis
            //                , normalIToflowD is parallel to x axis
            //                , normalIIToflowD is parallel to y axis
	    flowD_.x() = 0;
	    flowD_.y() = 0;
	    flowD_.z() = 1;

	    normalToflowDI_.x() = 1;
	    normalToflowDI_.y() = 0;
	    normalToflowDI_.z() = 0;

	    normalToflowDII_.x() = 0;
	    normalToflowDII_.y() = 1;
	    normalToflowDII_.z() = 0;
        }   

   	bDataL_[bDI].pressureDragForce = bDataL_[bDI].pressureForce & flowD_;
   	bDataL_[bDI].frictionDragForce = bDataL_[bDI].shearStressForce & flowD_;
        bDataL_[bDI].totalDragForce = bDataL_[bDI].pressureDragForce + bDataL_[bDI].frictionDragForce;

   	bDataL_[bDI].pressureLiftForceI = bDataL_[bDI].pressureForce & normalToflowDI_;
   	bDataL_[bDI].pressureLiftForceII = bDataL_[bDI].pressureForce & normalToflowDII_;
   	bDataL_[bDI].frictionLiftForceI = bDataL_[bDI].shearStressForce & normalToflowDI_;
   	bDataL_[bDI].frictionLiftForceII = bDataL_[bDI].shearStressForce & normalToflowDII_;

        bDataL_[bDI].totalLiftForceI = bDataL_[bDI].pressureLiftForceI + bDataL_[bDI].frictionLiftForceI;   
        bDataL_[bDI].totalLiftForceII = bDataL_[bDI].pressureLiftForceII + bDataL_[bDI].frictionLiftForceII;

   	bDataL_[bDI].velocity = (bDataL_[bDI].centre - bDataL_[bDI].centreOld)/(time-timeOld_+ROOTVSMALL);
   	bDataL_[bDI].accelaration = (bDataL_[bDI].velocity - bDataL_[bDI].velocityOld)/(time-timeOld_+ROOTVSMALL);
  	bDataL_[bDI].velocityOld = bDataL_[bDI].velocity;
   	bDataL_[bDI].centreOld = bDataL_[bDI].centre;
        timeOld_ = time;

        scalar sphereBubbleSurfaceArea = mathematical::pi * sqr( bDataL_[bDI].sphereBubbleDiameter);
        scalar magVelocityDiff = mag( magFlowVelocity_ * flowD_ - bDataL_[bDI].velocity);

   	bDataL_[bDI].dragCoeff = bDataL_[bDI].totalDragForce
                     /(0.5 * bDataL_[bDI].outerFluidDensity * (pow(magVelocityDiff,2.0)+ROOTVSMALL) * sphereBubbleSurfaceArea);

   	bDataL_[bDI].LiftCoeffI = bDataL_[bDI].totalLiftForceI
                     /(0.5 * bDataL_[bDI].outerFluidDensity * (pow(magVelocityDiff,2.0)+ROOTVSMALL) * sphereBubbleSurfaceArea);

   	bDataL_[bDI].LiftCoeffII = bDataL_[bDI].totalLiftForceII
                     /(0.5 * bDataL_[bDI].outerFluidDensity * (pow(magVelocityDiff,2.0)+ROOTVSMALL) * sphereBubbleSurfaceArea);

       // calculation of the droplet’s maximum extension in each coordinate direction, i.e
       // flowD_, normalToflowDI_, and normalToflowDII_.
       scalar maxDiameterComponent1 = 0.0; // it's in the direction of flowD_
       scalar maxDiameterComponent2 = 0.0; // it's in the direction of normalToflowDI_
       scalar maxDiameterComponent3 = 0.0; // it's in the direction of normalToflowDII_
       for(int i = 0; i < pL.size() ; i++)
       {
           for(int j = i+1; j < pL.size() ; j++)
           {
               vector frontPointToPoint = pL[i].currentPoint - pL[j].currentPoint;
               scalar diameterComponent1 = frontPointToPoint & flowD_; 
               scalar diameterComponent2 = frontPointToPoint & normalToflowDI_; 
               scalar diameterComponent3 = frontPointToPoint & normalToflowDII_; 

               if (  diameterComponent1 > maxDiameterComponent1)
               {
                   maxDiameterComponent1 = diameterComponent1;
               }
               if (  diameterComponent2 > maxDiameterComponent2)
               {
                   maxDiameterComponent2 = diameterComponent2;
               }
               if (  diameterComponent3 > maxDiameterComponent3)
               {
                   maxDiameterComponent3 = diameterComponent3;
               }
           }
       }
       bDataL_[bDI].diameter[0] = maxDiameterComponent1;
       bDataL_[bDI].diameter[1] = maxDiameterComponent2;
       bDataL_[bDI].diameter[2] = maxDiameterComponent3;
    }//1     
}

//note: func to print the Bubble Positions for all times.
template<class CloudType>
inline void Foam::FrontDynamic<CloudType>::writeBubblePositions()
{
    Info << "\n"  << "----------> Now writeBubblePositions is done ........ " << endl;

    fileName outfileName = "BubblePositionsPLT";

    // This is the directory keeping the bubble data. 
    fileName FTResult = "FTResult";

    fileName outfileNamePath=this->owner().db().time().rootPath()/this->owner().db().time().globalCaseName()/FTResult/outfileName;

    ofstream outFile(outfileNamePath.c_str(), ios::app);
    ostream& os = outFile;

    label bubbleNum = bDataL_.size();

    const scalar time = this->owner().db().time().value();
    const scalar deltaTime = this->owner().db().time().deltaTValue();

    if( (time-deltaTime) == timeTobeSteady_ )
    //if( firstTimeWriting_ == 1 )
    {  
        os  << " VARIABLES = \"time\" ";
        for(int bI = 0; bI < bubbleNum ; bI++)
        {//0
           OStringStream ossX;
           ossX << "bubble_" << bI << "_CentreX";
           fileName bCentreX = ossX.str();

           OStringStream ossY;
           ossY << "bubble_" << bI << "_CentreY";
           fileName bCentreY = ossY.str();

           OStringStream ossZ;
           ossZ << "bubble_" << bI << "_CentreZ";
           fileName bCentreZ = ossZ.str();

           os  << ",\"" << bCentreX << "\", \"" << bCentreY << "\", \"" <<bCentreZ << "\"";
        }//0
        os  << " \n";
    }

    os  << time << setw(15);  
    for(int bI = 0; bI < bubbleNum ; bI++)
    {//0
	os  << bDataL_[bI].centre.x() << setw(15)
	    << bDataL_[bI].centre.y() << setw(15)
	    << bDataL_[bI].centre.z() << token::SPACE;
    }//0
    os << " \n";

    outFile.close();
}

//note: func to print the Bubble Positions for all times.
template<class CloudType>
inline void Foam::FrontDynamic<CloudType>::writeBubbleVelocityAndAccelaration()
{
    Info << "\n"  << "----------> Now writeBubbleVelocityAndAccelaration is done ........ " << endl;

    fileName outfileName = "BubbleVelocityAndAccelarationPLT";

    // This is the directory keeping the bubble data. 
    fileName FTResult = "FTResult";

    fileName outfileNamePath=this->owner().db().time().rootPath()/this->owner().db().time().globalCaseName()/FTResult/outfileName;

    ofstream outFile(outfileNamePath.c_str(), ios::app);
    ostream& os = outFile;

    label bubbleNum = bDataL_.size();

    const scalar time = this->owner().db().time().value();
    const scalar deltaTime = this->owner().db().time().deltaTValue();

    if( (time-deltaTime) == timeTobeSteady_ )
    {  
        os  << " VARIABLES = \"time\" ";
        for(int bI = 0; bI < bubbleNum ; bI++)
        {//0

           OStringStream ossBV;
           ossBV << "bulkVelocity";
           fileName bBulkVelocity = ossBV.str();

           OStringStream ossUX;
           ossUX << "bubble_" << bI << "_CentreUX";
           fileName bCentreUX = ossUX.str();

           OStringStream ossUY;
           ossUY << "bubble_" << bI << "_CentreUY";
           fileName bCentreUY = ossUY.str();

           OStringStream ossUZ;
           ossUZ << "bubble_" << bI << "_CentreUZ";
           fileName bCentreUZ = ossUZ.str();

           OStringStream ossAX;
           ossAX << "bubble_" << bI << "_CentreAX";
           fileName bCentreAX = ossAX.str();

           OStringStream ossAY;
           ossAY << "bubble_" << bI << "_CentreAY";
           fileName bCentreAY = ossAY.str();

           OStringStream ossAZ;
           ossAZ << "bubble_" << bI << "_CentreAZ";
           fileName bCentreAZ = ossAZ.str();

           os  << ",\"" << bBulkVelocity << "\", \"" <<  bCentreUX << "\", \"" << bCentreUY << "\", \"" <<bCentreUZ << "\", "
               << "\"" << bCentreAX << "\", \"" << bCentreAY << "\", \"" <<bCentreAZ << "\"";
        }//0
        os  << " \n";
    }

    os  << time << setw(15);  
    for(int bI = 0; bI < bubbleNum ; bI++)
    {//0
	os  << magFlowVelocity_ << setw(15)
	    << bDataL_[bI].velocity.x() << setw(15)
	    << bDataL_[bI].velocity.y() << setw(15)
	    << bDataL_[bI].velocity.z() << setw(15)
	    << bDataL_[bI].accelaration.x() << setw(15)
	    << bDataL_[bI].accelaration.y() << setw(15)
	    << bDataL_[bI].accelaration.z() << token::SPACE;
    }//0
    os << " \n";

    outFile.close();
}

//note: func to print the Bubble Positions for all times.
template<class CloudType>
inline void Foam::FrontDynamic<CloudType>::writeBubbleSphericity()
{
    Info << "\n"  << "----------> Now writeBubbleSphericity is done ........ " << endl;

    fileName outfileName = "BubbleSphericityPLT";

    // This is the directory keeping the bubble data. 
    fileName FTResult = "FTResult";

    fileName outfileNamePath=this->owner().db().time().rootPath()/this->owner().db().time().globalCaseName()/FTResult/outfileName;

    ofstream outFile(outfileNamePath.c_str(), ios::app);
    ostream& os = outFile;

    label bubbleNum = bDataL_.size();

    const scalar time = this->owner().db().time().value();
    const scalar deltaTime = this->owner().db().time().deltaTValue();

    if( (time-deltaTime) == timeTobeSteady_ )
    {  
        os  << " VARIABLES = \"time\" ";
        for(int bI = 0; bI < bubbleNum ; bI++)
        {//0
           OStringStream ossSphericity;
           ossSphericity << "bubble_" << bI << "_Sphericity";
           fileName bSphericity = ossSphericity.str();

           OStringStream ossSurfaceArea;
           ossSurfaceArea << "bubble_" << bI << "_SurfaceArea";
           fileName bSurfaceArea = ossSurfaceArea.str();

           OStringStream ossVolume;
           ossVolume << "bubble_" << bI << "_Volume";
           fileName bVolume = ossVolume.str();

           os  << ",\"" << bSphericity << "\", \"" << bSurfaceArea << "\", \"" << bVolume << "\"";
        }//0
        os  << " \n";
    }

    os  << time << setw(15);  
    for(int bI = 0; bI < bubbleNum ; bI++)
    {//0
	os  << bDataL_[bI].sphericity << setw(15)
	    << bDataL_[bI].surfaceArea << setw(15)
	    << bDataL_[bI].volume << token::SPACE;
    }//0
    os << " \n";

    outFile.close();
}

//note: func to print the Bubble Positions for all times.
template<class CloudType>
inline void Foam::FrontDynamic<CloudType>::writeBubbleDragCoefs()
{
    Info << "\n"  << "----------> Now writeBubbleDragCoefs is done ........ " << endl;

    fileName outfileName = "BubbleDragCoefsPLT";

    // This is the directory keeping the bubble data. 
    fileName FTResult = "FTResult";

    fileName outfileNamePath=this->owner().db().time().rootPath()/this->owner().db().time().globalCaseName()/FTResult/outfileName;

    ofstream outFile(outfileNamePath.c_str(), ios::app);
    ostream& os = outFile;

    label bubbleNum = bDataL_.size();

    const scalar time = this->owner().db().time().value();
    const scalar deltaTime = this->owner().db().time().deltaTValue();

    if( (time-deltaTime) == timeTobeSteady_ )
    {  
        os  << " VARIABLES = \"time\" ";
        for(int bI = 0; bI < bubbleNum ; bI++)
        {//0
           OStringStream ossPDC;
           ossPDC << "bubble_" << bI << "_PressureDragCoef";
           fileName bPressureDragCoef = ossPDC.str();

           OStringStream ossFDC;
           ossFDC << "bubble_" << bI << "_FrictionDragCoef";
           fileName bFrictionDragCoef = ossFDC.str();

           OStringStream ossTDC;
           ossTDC << "bubble_" << bI << "_TotalDragCoef";
           fileName bTotalDragCoef = ossTDC.str();

           os  << ",\"" << bPressureDragCoef << "\", \"" << bFrictionDragCoef << "\", \"" << bTotalDragCoef << "\"";
        }//0
        os  << " \n";
    }

    os  << time << setw(15);  
    for(int bI = 0; bI < bubbleNum ; bI++)
    {//0
        scalar denaminator = bDataL_[bI].totalDragForce/bDataL_[bI].dragCoeff + ROOTVSMALL;
	os  << bDataL_[bI].pressureDragForce/denaminator << setw(15)
	    << bDataL_[bI].frictionDragForce/denaminator << setw(15)
	    << bDataL_[bI].dragCoeff << token::SPACE;
    }//0
    os << " \n";

    outFile.close();
}

//note: func to print the Bubble Positions for all times.
template<class CloudType>
inline void Foam::FrontDynamic<CloudType>::writeBubbleLiftCoefs()
{
    Info << "\n"  << "----------> Now writeBubbleLiftCoefs is done ........ " << endl;

    fileName outfileName = "BubbleLiftCoefsPLT";

    // This is the directory keeping the bubble data. 
    fileName FTResult = "FTResult";

    fileName outfileNamePath=this->owner().db().time().rootPath()/this->owner().db().time().globalCaseName()/FTResult/outfileName;

    ofstream outFile(outfileNamePath.c_str(), ios::app);
    ostream& os = outFile;

    label bubbleNum = bDataL_.size();

    const scalar time = this->owner().db().time().value();
    const scalar deltaTime = this->owner().db().time().deltaTValue();

    if( (time-deltaTime) == timeTobeSteady_ )
    {  
        os  << " VARIABLES = \"time\" ";
        for(int bI = 0; bI < bubbleNum ; bI++)
        {//0
           OStringStream ossPLCI;
           ossPLCI << "bubble_" << bI << "_PressureLiftCoefI";
           fileName bPressureLiftCoefI = ossPLCI.str();

           OStringStream ossPLCII;
           ossPLCII << "bubble_" << bI << "_PressureLiftCoefII";
           fileName bPressureLiftCoefII = ossPLCII.str();

           OStringStream ossFLCI;
           ossFLCI << "bubble_" << bI << "_FrictionLiftCoefI";
           fileName bFrictionLiftCoefI = ossFLCI.str();

           OStringStream ossFLCII;
           ossFLCII << "bubble_" << bI << "_FrictionLiftCoefII";
           fileName bFrictionLiftCoefII = ossFLCII.str();

           OStringStream ossTLCI;
           ossTLCI << "bubble_" << bI << "_TotalLiftCoefI";
           fileName bTotalLiftCoefI = ossTLCI.str();

           OStringStream ossTLCII;
           ossTLCII << "bubble_" << bI << "_TotalLiftCoefII";
           fileName bTotalLiftCoefII = ossTLCII.str();

           os  << ",\"" << bPressureLiftCoefI << "\", \"" << bPressureLiftCoefII << "\", \""
                        << bFrictionLiftCoefI << "\", \"" << bFrictionLiftCoefII << "\", \"" 
                        << bTotalLiftCoefI << "\", \"" << bTotalLiftCoefII << "\"";
        }//0
        os  << " \n";
    }

    os  << time << setw(15);  
    for(int bI = 0; bI < bubbleNum ; bI++)
    {//0
        scalar denaminator = bDataL_[bI].totalDragForce/bDataL_[bI].dragCoeff + ROOTVSMALL;
	os  << bDataL_[bI].pressureLiftForceI/denaminator << setw(15)
	    << bDataL_[bI].pressureLiftForceII/denaminator << setw(15)
	    << bDataL_[bI].frictionLiftForceI/denaminator << setw(15)
	    << bDataL_[bI].frictionLiftForceII/denaminator << setw(15)
	    << bDataL_[bI].LiftCoeffI << setw(15)
	    << bDataL_[bI].LiftCoeffII << token::SPACE;
    }//0
    os << " \n";

    outFile.close();
}

//note: func to print the Bubble Diameter Components.
template<class CloudType>
inline void Foam::FrontDynamic<CloudType>::writeBubbleDiameterComponents()
{
    Info << "\n"  << "----------> Now writeBubbleDiameterComponents is done ........ " << endl;

    fileName outfileName = "BubbleDiameterComponentsPLT";

    // This is the directory keeping the bubble data. 
    fileName FTResult = "FTResult";

    fileName outfileNamePath=this->owner().db().time().rootPath()/this->owner().db().time().globalCaseName()/FTResult/outfileName;

    ofstream outFile(outfileNamePath.c_str(), ios::app);
    ostream& os = outFile;

    label bubbleNum = bDataL_.size();

    const scalar time = this->owner().db().time().value();
    const scalar deltaTime = this->owner().db().time().deltaTValue();

    if( (time-deltaTime) == timeTobeSteady_ )
    {  
        os  << " VARIABLES = \"time\" ";
        for(int bI = 0; bI < bubbleNum ; bI++)
        {//0
           OStringStream ossDI;
           ossDI << "bubble_" << bI << "_DI";
           fileName bDI = ossDI.str();

           OStringStream ossDII;
           ossDII << "bubble_" << bI << "_DII";
           fileName bDII = ossDII.str();

           OStringStream ossDIII;
           ossDIII << "bubble_" << bI << "_DIII";
           fileName bDIII = ossDIII.str();

           os  << ",\"" << bDI << "\", \"" << bDII << "\", \"" << bDIII << "\"";
        }//0
        os  << " \n";
    }

    os  << time << setw(15);  
    for(int bI = 0; bI < bubbleNum ; bI++)
    {//0
	os  << bDataL_[bI].diameter[0] << setw(15)
	    << bDataL_[bI].diameter[1] << setw(15)
	    << bDataL_[bI].diameter[2] << token::SPACE;
    }//0
    os << " \n";

    outFile.close();
}

//note: func to calculate non-dimensional numbers.
template<class CloudType>
inline void Foam::FrontDynamic<CloudType>::calcBubbleNonDimensionalNumbers()
{
    scalar gravity = mag (this->owner().g().value() );
    scalar basedDensity = bDataL_[0].outerFluidDensity;
    scalar basedViscosity = bDataL_[0].outerFluidViscosity;

    forAll(bDataL_,bDI)
    {//1
        scalar db = bDataL_[bDI].sphereBubbleDiameter;
        bDataL_[bDI].surfaceTensionCoeff = basedDensity*pow(db, 2.0)*gravity/bDataL_[bDI].EotvosNumber;
        scalar sigma =bDataL_[bDI].surfaceTensionCoeff;
        bDataL_[bDI].MortonNumber = gravity*pow(basedViscosity, 4)/basedDensity/pow(sigma, 3);
    }
}

//note: func to do all postprocessing jobs.
template<class CloudType>
inline void Foam::FrontDynamic<CloudType>::outputTimePostProcessing()
{
    //if( this->owner().db().time().outputTime() && Pstream::master() )
    {
        Info << "\n----------> At the time: " << this->owner().db().time().value()
             <<" , Now all outputTimePostProcessing jobs is doing ....... " << endl;
        bubblePostProcessing();
        writeBubblePositions();
        writeBubbleVelocityAndAccelaration();
        writeBubbleSphericity();
        writeBubbleDragCoefs();
        writeBubbleLiftCoefs();
        writeBubbleDiameterComponents();
        printFrontsForAllTimes();
        writeDataOfTheFrontsAtTheOutputTime();
        writeDataOfTheFrontsAtTheOutputTimeIntimeName();
    }
}

//note: func to calculate flow direction field.
template<class CloudType>
inline void Foam::FrontDynamic<CloudType>::flowDirectionField()
{
    const fvMesh& mesh = this->owner().mesh();
    const pointField& ctrs = mesh.cellCentres();

    // calculating the flow direction:
    if( periodicalOption_ == "rotationalPeriodicalInCurvedDuct" || periodicalOption_ == "combinedRotationalPeriodicalInCurvedDuct" )
    {
        forAll(ctrs, cellI)
        {
	    rCentre_.z() =  ctrs[cellI].z();
	    vector rStar = ctrs[cellI] - rCentre_;

	    flowD_.x() = -rStar.y();
	    flowD_.y() =  rStar.x();
	    flowD_.z() =  0;
	    flowD_ = flowD_/mag(flowD_);
            this->owner().dPressureFromFT()[cellI] = flowD_ * gradPAtflowDirection_;
            flowDList_.append(flowD_);
        }   
    } 
    if( periodicalOption_ == "translationalPeriodicalInCurvedDuct" || periodicalOption_ == "combinedTranslationalPeriodicalInCurvedDuct")
    {
        forAll(ctrs, cellI)
        {
	    flowD_.x() = 0;
	    flowD_.y() = 0;
	    flowD_.z() = 1;
            this->owner().dPressureFromFT()[cellI] = flowD_ * gradPAtflowDirection_;
            flowDList_.append(flowD_);
        }   
    }
    else if( periodicalOption_ == "straightDuctPeriodical" || periodicalOption_ == "none" )
    {
        // for none option, flowD is parallel to z axis
        //                , normalIToflowD is parallel to x axis
        //                , normalIIToflowD is parallel to y axis
        forAll(ctrs, cellI)
        {
	    flowD_.x() = 0;
	    flowD_.y() = 0;
	    flowD_.z() = 1;
            this->owner().dPressureFromFT()[cellI] = flowD_ * gradPAtflowDirection_;
            flowDList_.append(flowD_);
        }  
    }
}

//note: func to print all important runTime data.
template<class CloudType>
inline void Foam::FrontDynamic<CloudType>::printingAllImportantRunTimeData()
{
    const fvMesh& mesh = this->owner().mesh();
    const pointField& ctrs = mesh.cellCentres();

    // Integrate flow variables over cell set
    scalar magUbarAve = 0.0;
    scalar rhoVol = 0.0;//Jafari added
    const scalarField& cv = mesh.V();
    forAll(ctrs, i)
    {
        label cellI = i;//cells_[i];
        scalar volCell = cv[cellI];
        rhoVol += densityIndicator_.internalField()[cellI] * volCell;//Jafari added
        magUbarAve += (flowDList_[cellI] & this->owner().U().internalField()[cellI])
                      *densityIndicator_.internalField()[cellI]*volCell;//Jafari edited
    }

    // Collect across all processors
    reduce(magUbarAve, sumOp<scalar>());
    reduce(rhoVol, sumOp<scalar>());//Jafari added

    // Volume averages
    magUbarAve /= rhoVol;//Jafari edited

    magFlowVelocity_ = magUbarAve;

    Info << "\n----------> At the time: " << this->owner().db().time().value()
         << " due to the pressure gradient value of " << gradPAtflowDirection_ 
         << " the average velocity at the flow direction is " << magUbarAve << endl;

    scalar basedDensity = bDataL_[0].outerFluidDensity;
    scalar basedViscosity = bDataL_[0].outerFluidViscosity;
    ReDh_ = basedDensity * magFlowVelocity_* Dh_/basedViscosity;
    Info << "\n----------> At the time: " << this->owner().db().time().value()
         << " the reynold number, ReDh_ is calculted as " << ReDh_ << endl;
}

//note: func to initialize the CPTList_.
template<class CloudType>
inline void Foam::FrontDynamic<CloudType>::initialAllocatingTheCPTList()
{
    //Info << "\n----------> At the time: " << this->owner().db().time().value()
    //     << " the CPTList_ is initially allocated. " << endl;

    const fvMesh& mesh = this->owner().mesh();
    const pointField& ctrs = mesh.cellCentres();

    CPTList_.clear();

    CPTList_.setSize(ctrs.size());

    forAll(CPTList_, i)
    {
        CPTList_[i].inducedDensity.setSize(bDataL_.size());
        CPTList_[i].inducedViscosity.setSize(bDataL_.size());
    }
}

//note: func to initialize the CPTList_.
template<class CloudType>
inline void Foam::FrontDynamic<CloudType>::setLengthScalesNearTheFront()
{
    //Info << "\n----------> At the time: " << this->owner().db().time().value()
    //     << " the CPTList_ is initially allocated. " << endl;

    lengthScaleNearTheFront_.setSize(bDataL_.size());

    forAll(bDataL_,bDI)
    {//1
        lengthScaleNearTheFront_[bDI] = h_;
    }//1
}
