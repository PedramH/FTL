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

// * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * //

inline Foam::torqueSuSp::torqueSuSp()
{}


inline Foam::torqueSuSp::torqueSuSp
(
    const Tuple2<vector, scalar>& fs
)
:
    Tuple2<vector, scalar>(fs)
{}


inline Foam::torqueSuSp::torqueSuSp(const vector& Su, const scalar Sp)
{
    first() = Su;
    second() = Sp;
}


inline Foam::torqueSuSp::torqueSuSp(Istream& is)
:
    Tuple2<vector, scalar>(is)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::vector& Foam::torqueSuSp::Su() const
{
    return first();
}


inline Foam::scalar Foam::torqueSuSp::Sp() const
{
    return second();
}


inline Foam::vector& Foam::torqueSuSp::Su()
{
    return first();
}


inline Foam::scalar& Foam::torqueSuSp::Sp()
{
    return second();
}


// * * * * * * * * * * * * * * * Operators * * * * * * * * * * * * * * * * * //

inline void Foam::torqueSuSp::operator=(const torqueSuSp& susp)
{
    first() = susp.first();
    second() = susp.second();
}


inline void Foam::torqueSuSp::operator+=(const torqueSuSp& susp)
{
    first() += susp.first();
    second() += susp.second();
}


inline void Foam::torqueSuSp::operator-=(const torqueSuSp& susp)
{
    first() -= susp.first();
    second() -= susp.second();
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

inline Foam::torqueSuSp Foam::operator+
(
    const torqueSuSp& susp1,
    const torqueSuSp& susp2
)
{
    return torqueSuSp
    (
        susp1.first() + susp2.first(),
        susp1.second() + susp2.second()
    );
}


inline Foam::torqueSuSp Foam::operator*
(
    const scalar s,
    const torqueSuSp& susp
)
{
    return torqueSuSp(susp.first()*s, susp.second()*s);
}


// ************************************************************************* //
