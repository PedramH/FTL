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

\*---------------------------------------------------------------------------*/

#include "mykinematicParcelInjectionData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::mykinematicParcelInjectionData::mykinematicParcelInjectionData(Istream& is)
:
    kinematicParcelInjectionData(is)
{
    is.check("reading deltaD");
    is >> deltaD_;

	is.check("reading uiuj");
    is >> uiuj_;

	is.check("reading kf");
    is >> kf_;

    is.check("mykinematicParcelInjectionData(Istream& is)");
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const mykinematicParcelInjectionData& data
)
{
    os << static_cast<const kinematicParcelInjectionData&>(data);

    os << data.deltaD_;

	os << data.uiuj_;

	os << data.kf_;

    return os;
}


Foam::Istream& Foam::operator>>(Istream& is, mykinematicParcelInjectionData& data)
{
    is >> static_cast<kinematicParcelInjectionData&>(data);

    is.check("reading deltaD");
    is >> data.deltaD_;

	is.check("reading uiuj");
    is >> data.uiuj_;

	is.check("reading kf");
    is >> data.kf_;

    is.check("operator(Istream&, mykinematicParcelInjectionData&)");

    return is;
}


// ************************************************************************* //
