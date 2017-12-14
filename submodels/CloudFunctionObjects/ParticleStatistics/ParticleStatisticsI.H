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

template<class CloudType>
inline const Foam::label&
Foam::ParticleStatistics<CloudType>::gNum() const
{
   return gNum_;
}

template<class CloudType>
inline const Foam::Switch&
Foam::ParticleStatistics<CloudType>::UMean() const
{
   return UMean_;
}

template<class CloudType>
inline const Foam::Switch&
Foam::ParticleStatistics<CloudType>::uPrime2Mean() const
{
   return uPrime2Mean_;
}

template<class CloudType>
inline const Foam::Switch&
Foam::ParticleStatistics<CloudType>::writed10() const
{
   return writed10_;
}

template<class CloudType>
inline const Foam::Switch&
Foam::ParticleStatistics<CloudType>::writed32() const
{
   return writed32_;
}

template<class CloudType>
inline const Foam::Switch&
Foam::ParticleStatistics<CloudType>::resetFields() const
{
   return resetFields_;
}

template<class CloudType>
Foam::label Foam::ParticleStatistics<CloudType>::timeStart() const
{
   return timeStart_;
}

template<class CloudType>
Foam::label Foam::ParticleStatistics<CloudType>::timeEnd() const
{
   return timeEnd_;
}



// ************************************************************************* //
