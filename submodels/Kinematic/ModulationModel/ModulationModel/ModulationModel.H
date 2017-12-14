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
    Foam::ModulationModel

Description

\*---------------------------------------------------------------------------*/

#ifndef ModulationModel_H
#define ModulationModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "CloudSubModelBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class ModulationModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class ModulationModel
:
    public CloudSubModelBase<CloudType>
{
public:

    //- Runtime type information
    TypeName("modulationModel");


    // Declare runtime constructor selection table

    declareRunTimeSelectionTable
    (
        autoPtr,
        ModulationModel,
        dictionary,
        (
            const dictionary& dict,
            CloudType& owner
        ),
        (dict, owner)
    );


    // Constructors

        //- Construct null from owner
        ModulationModel(CloudType& owner);

        //- Construct from components
        ModulationModel
        (
            const dictionary& dict,
            CloudType& owner,
            const word& type
        );

        //- Construct copy
        ModulationModel(ModulationModel<CloudType>& dm);

        //- Construct and return a clone
        virtual autoPtr<ModulationModel<CloudType> > clone()
        {
            return autoPtr<ModulationModel<CloudType> >
            (
                new ModulationModel<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~ModulationModel();


    //- Selector
    static autoPtr<ModulationModel<CloudType> > New
    (
        const dictionary& dict,
        CloudType& owner
    );


    // Member Functions

        //- Update (disperse particles)
        virtual symmTensor update
		(
			const vector UcMean,
			const vector Uc,
			const vector Up0,
			const vector Upn,
			const vector Su
			//const typename CloudType::parcelType& p
		);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeModulationModel(CloudType)                                        \
                                                                              \
    typedef CloudType::kinematicCloudType kinematicCloudType;                 \
    defineTemplateTypeNameAndDebug(ModulationModel<kinematicCloudType>, 0);   \
    defineTemplateRunTimeSelectionTable                                       \
    (                                                                         \
        ModulationModel<kinematicCloudType>,                                  \
        dictionary                                                            \
    );


#define makeModulationModelType(SS, CloudType)                                \
                                                                              \
    typedef CloudType::kinematicCloudType kinematicCloudType;                 \
    defineNamedTemplateTypeNameAndDebug(SS<kinematicCloudType>, 0);           \
                                                                              \
    ModulationModel<kinematicCloudType>::                                     \
        adddictionaryConstructorToTable<SS<kinematicCloudType> >              \
            add##SS##CloudType##kinematicCloudType##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "ModulationModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
