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

#include "ElectrokineticParcel.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class ParcelType>
Foam::string Foam::ElectrokineticParcel<ParcelType>::propertyList_ =
    Foam::ElectrokineticParcel<ParcelType>::propertyList();


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::ElectrokineticParcel<ParcelType>::ElectrokineticParcel
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields
)
:
    ParcelType(mesh, is, readFields),
    z_(0),
	Ec_(vector::zero)
{
    if (readFields)
    {
        if (is.format() == IOstream::ASCII)
        {
            is >> z_;
        }
        else
        {
            is.read
            (
                reinterpret_cast<char*>(&z_),
              + sizeof(z_)
            );
        }
    }

    // Check state of Istream
    is.check
    (
        "ElectrokineticParcel::ElectrokineticParcel(const polyMesh&, Istream&, bool)"
    );
}


template<class ParcelType>
template<class CloudType>
void Foam::ElectrokineticParcel<ParcelType>::readFields(CloudType& c)
{
    if (!c.size())
    {
        return;
    }

    ParcelType::readFields(c);

    IOField<scalar> z(c.fieldIOobject("z", IOobject::MUST_READ));
    c.checkFieldIOobject(c, z);

    label i = 0;
    forAllIter(typename Cloud<ElectrokineticParcel<ParcelType> >, c, iter)
    {
        ElectrokineticParcel<ParcelType>& p = iter();

        p.z_ = z[i];
        i++;
    }
}


template<class ParcelType>
template<class CloudType>
void Foam::ElectrokineticParcel<ParcelType>::writeFields(const CloudType& c)
{
    ParcelType::writeFields(c);

    label np = c.size();

   	IOField<scalar> z(c.fieldIOobject("z", IOobject::NO_READ), np);
    
    label i = 0;
    forAllConstIter(typename Cloud<ElectrokineticParcel<ParcelType> >, c, iter)
    {
        const ElectrokineticParcel<ParcelType>& p = iter();

        z[i] = p.z_;
        i++;
    }

    z.write();
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class ParcelType>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const ElectrokineticParcel<ParcelType>& p
)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << static_cast<const ParcelType&>(p)
            << token::SPACE << p.z();
    }
    else
    {
        os  << static_cast<const ParcelType&>(p);
        os.write
        (
            reinterpret_cast<const char*>(&p.z_),
            sizeof(p.z()) 
        );
    }

    // Check state of Ostream
    os.check
    (
        "Ostream& operator<<(Ostream&, const ElectrokineticParcel<ParcelType>&)"
    );

    return os;
}


// ************************************************************************* //
