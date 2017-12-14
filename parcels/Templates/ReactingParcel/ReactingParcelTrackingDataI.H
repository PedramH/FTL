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

template<class ParcelType>
template<class CloudType>
inline Foam::ReactingParcel<ParcelType>::TrackingData<CloudType>::TrackingData
(
    CloudType& cloud,
    trackPart part,
	fieldUpdate updates //added1
)
:
    ParcelType::template TrackingData<CloudType>(cloud, part, updates), //added1
    pInterp_
    (
        interpolation<scalar>::New
        (
            cloud.solution().interpolationSchemes(),
            cloud.p()
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ParcelType>
template<class CloudType>
Foam::ReactingParcel<ParcelType>::TrackingData<CloudType>::~TrackingData()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParcelType> //added1
template<class CloudType>
inline void
Foam::ReactingParcel<ParcelType>::TrackingData<CloudType>::updateParticleStatistics
(
    //CloudType& cloud
)
{
	ParcelType::template TrackingData<CloudType>::updateParticleStatistics(); //(cloud);
	//CloudType& cloud = this->cloud();
}


template<class ParcelType>
template<class CloudType>
inline const Foam::interpolation<Foam::scalar>&
Foam::ReactingParcel<ParcelType>::TrackingData<CloudType>::pInterp() const
{
    return pInterp_();
}


// ************************************************************************* //
