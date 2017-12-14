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
    Foam::NoFrontAnalyzing

Description
    Dummy collision model for 'none'

\*---------------------------------------------------------------------------*/

#ifndef NoFrontAnalyzing_H
#define NoFrontAnalyzing_H

#include "FrontAnalyzingModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                    Class NoFrontAnalyzing Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class NoFrontAnalyzing
:
    public FrontAnalyzingModel<CloudType>
{
protected:

    // Protected Member Functions

        //- Update the model
        virtual void collide(const scalar dt);


public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct from dictionary
        NoFrontAnalyzing(const dictionary& dict, CloudType& owner);

        //- Construct copy
        NoFrontAnalyzing(NoFrontAnalyzing<CloudType>& ft);

        //- Construct and return a clone
        virtual autoPtr<FrontAnalyzingModel<CloudType> > clone() //const
        {
            return autoPtr<FrontAnalyzingModel<CloudType> >
            (
                new NoFrontAnalyzing<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~NoFrontAnalyzing();


    // Member Functions

        //- Flag to indicate whether model activates collision model
        virtual bool active() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "NoFrontAnalyzing.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
