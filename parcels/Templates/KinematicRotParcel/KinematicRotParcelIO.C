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

#include "KinematicRotParcel.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class ParcelType>
Foam::string Foam::KinematicRotParcel<ParcelType>::propertyList_ =
    Foam::KinematicRotParcel<ParcelType>::propertyList();


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::KinematicRotParcel<ParcelType>::KinematicRotParcel
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields
)
:
    ParcelType(mesh, is, readFields),
    Omega_(vector::zero),
	Omegac_(vector::zero)
{
    if (readFields)
    {
        if (is.format() == IOstream::ASCII)
        {
            is >> Omega_;
        }
        else
        {
            is.read
            (
                reinterpret_cast<char*>(&Omega_),
              + sizeof(Omega_)
            );
        }
    }

    // Check state of Istream
    is.check
    (
        "KinematicRotParcel::KinematicRotParcel(const polyMesh&, Istream&, bool)"
    );
}


template<class ParcelType>
template<class CloudType>
void Foam::KinematicRotParcel<ParcelType>::readFields(CloudType& c)
{
    if (!c.size())
    {
        return;
    }

    ParcelType::readFields(c);

    IOField<vector> Omega(c.fieldIOobject("Omega", IOobject::MUST_READ));
    c.checkFieldIOobject(c, Omega);

    label i = 0;
    forAllIter(typename Cloud<KinematicRotParcel<ParcelType> >, c, iter)
    {
        KinematicRotParcel<ParcelType>& p = iter();

        p.Omega_ = Omega[i];
        i++;
    }
}


template<class ParcelType>
template<class CloudType>
void Foam::KinematicRotParcel<ParcelType>::writeFields(const CloudType& c)
{
    ParcelType::writeFields(c);

    label np = c.size();

   	IOField<vector> Omega(c.fieldIOobject("Omega", IOobject::NO_READ), np);
    
    label i = 0;
    forAllConstIter(typename Cloud<KinematicRotParcel<ParcelType> >, c, iter)
    {
        const KinematicRotParcel<ParcelType>& p = iter();

        Omega[i] = p.Omega_;
        i++;
    }

    Omega.write();
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class ParcelType>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const KinematicRotParcel<ParcelType>& p
)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << static_cast<const ParcelType&>(p)
            << token::SPACE << p.Omega();
    }
    else
    {
        os  << static_cast<const ParcelType&>(p);
        os.write
        (
            reinterpret_cast<const char*>(&p.Omega_),
            sizeof(p.Omega()) 
        );
    }

    // Check state of Ostream
    os.check
    (
        "Ostream& operator<<(Ostream&, const KinematicRotParcel<ParcelType>&)"
    );

    return os;
}


// ************************************************************************* //
