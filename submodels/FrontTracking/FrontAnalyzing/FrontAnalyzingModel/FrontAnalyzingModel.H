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
    Foam::FrontAnalyzingModel

Description
    Templated stochastic collision model class

SourceFiles
    FrontAnalyzingModel.C
    FrontAnalyzingModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef FrontAnalyzingModel_H
#define FrontAnalyzingModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "CloudSubModelBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class FrontAnalyzingModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class FrontAnalyzingModel
:
    public CloudSubModelBase<CloudType>
{
protected:

    //- Main collision routine
    virtual void dynamic();


public:

    //- Runtime type information
    TypeName("frontAnalyzingModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        FrontAnalyzingModel,
        dictionary,
        (
            const dictionary& dict,
            CloudType& owner
        ),
        (dict, owner)
    );


    // Constructors

        //- Construct null from owner
        FrontAnalyzingModel(CloudType& owner);

        //- Construct from dictionary
        FrontAnalyzingModel
        (
            const dictionary& dict,
            CloudType& owner,
            const word& type
        );

        //- Construct copy
        FrontAnalyzingModel(FrontAnalyzingModel<CloudType>& ft);

        //- Construct and return a clone
        virtual autoPtr<FrontAnalyzingModel<CloudType> > clone() //const
        {
            return autoPtr<FrontAnalyzingModel<CloudType> >
            (
                new FrontAnalyzingModel<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~FrontAnalyzingModel();


    //- Selector
    static autoPtr<FrontAnalyzingModel<CloudType> > New
    (
        const dictionary& dict,
        CloudType& owner
    );


    // Member Functions

        //- Update the model
        void update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#define makeFrontAnalyzingModel(CloudType)                                    \
                                                                              \
    typedef CloudType::frontTrackingCloudType frontTrackingCloudType;         \
    defineNamedTemplateTypeNameAndDebug                                       \
    (                                                                         \
        FrontAnalyzingModel<frontTrackingCloudType>,                          \
        0                                                                     \
    );                                                                        \
    defineTemplateRunTimeSelectionTable                                       \
    (                                                                         \
        FrontAnalyzingModel<frontTrackingCloudType>,                          \
        dictionary                                                            \
    );


#define makeFrontAnalyzingModelType(SS, CloudType)                            \
                                                                              \
    typedef CloudType::frontTrackingCloudType frontTrackingCloudType;         \
    defineNamedTemplateTypeNameAndDebug(SS<frontTrackingCloudType>, 0);       \
                                                                              \
    FrontAnalyzingModel<frontTrackingCloudType>::                             \
        adddictionaryConstructorToTable<SS<frontTrackingCloudType> >          \
            add##SS##CloudType##frontTrackingCloudType##ConstructorToTable_;



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "FrontAnalyzingModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
