/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::harmonic0

Group
    grpFvSurfaceInterpolationSchemes

Description
    harmonic0-mean differencing scheme class.

    This scheme interpolates 1/field using a scheme specified at run-time
    and return the reciprocal of the interpolate.

SourceFiles
    harmonic0.C

\*---------------------------------------------------------------------------*/

#ifndef harmonic0_H
#define harmonic0_H

#include "surfaceInterpolationScheme.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "reverseLinear.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class harmonic0 Declaration
\*---------------------------------------------------------------------------*/

class harmonic0
:
    public surfaceInterpolationScheme<scalar>
{
    // Private Member Functions

        //- No copy assignment
        void operator=(const harmonic0&) = delete;


public:

    //- Runtime type information
    TypeName("harmonic0");


    // Constructors

        //- Construct from mesh
        harmonic0(const fvMesh& mesh)
        :
            surfaceInterpolationScheme<scalar>(mesh)
        {}

        //- Construct from Istream.
        //  The name of the flux field is read from the Istream and looked-up
        //  from the mesh objectRegistry
        harmonic0
        (
            const fvMesh& mesh,
            Istream& is
        )
        :
            surfaceInterpolationScheme<scalar>(mesh)
        {}

        //- Construct from faceFlux and Istream
        harmonic0
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& is
        )
        :
            surfaceInterpolationScheme<scalar>(mesh)
        {}


    // Member Functions

        //- Return the interpolation weighting factors
        virtual tmp<surfaceScalarField> weights
        (
            const GeometricField<scalar, fvPatchField, volMesh>&
        ) const
        {
            NotImplemented;

            return tmp<surfaceScalarField>(nullptr);
        }

        //- Return the face-interpolate of the given cell field
        virtual tmp<GeometricField<scalar, fvsPatchField, surfaceMesh>>
        interpolate
        (
            const GeometricField<scalar, fvPatchField, volMesh>& vf
        ) const
        {
            return 1.0/(reverseLinear<scalar>(vf.mesh()).interpolate(1.0/(vf+dimensionedScalar("small", vf.dimensions(), VSMALL))));
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

