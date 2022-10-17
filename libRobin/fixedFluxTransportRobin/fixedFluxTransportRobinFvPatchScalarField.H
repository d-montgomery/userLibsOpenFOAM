/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2019
     \\/     M anipulation  | Matteo Icardi, Federico Municchi
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM.

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
    Foam::fixedFluxTransportRobinFvPatchScalarField

Description
    This boundary condition provide a Robin BC with an additional 
    flux read from a surfaceVectorField. Returns zero coefficients if 
    flux is zero.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type  fixedFluxTransportRobin;

        phi   phiv;

        RobinK  uniform 0;

        RobinD  uniform 1;

        RobinF  uniform 0;

        value $internalField;
    }
    \endverbatim


SourceFiles
    fixedFluxTransportRobinFvPatchScalarField.C

Authors:
    Federico Municchi (2022)
\*---------------------------------------------------------------------------*/

#ifndef fixedFluxTransportRobinFvPatchScalarField_H
#define fixedFluxTransportRobinFvPatchScalarField_H

#include "ghostRobinFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class fixedFluxTransportRobinFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class fixedFluxTransportRobinFvPatchScalarField
:
    public ghostRobinFvPatchScalarField
{
private:

    //- Flux variable
    word phiName_;

    //- use rho in diffusion coefficient
    bool useRho_;

    //- Effective RobinK_
    // scalarField RobinKeff_;

public:

    //- Runtime type information
    TypeName("fixedFluxTransportRobin");


    // Constructors

        //- Construct from patch and internal field
        fixedFluxTransportRobinFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fixedFluxTransportRobinFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given fixedFluxTransportRobinFvPatchScalarField
        //  onto a new patch
        fixedFluxTransportRobinFvPatchScalarField
        (
            const fixedFluxTransportRobinFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        fixedFluxTransportRobinFvPatchScalarField
        (
            const fixedFluxTransportRobinFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new fixedFluxTransportRobinFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        fixedFluxTransportRobinFvPatchScalarField
        (
            const fixedFluxTransportRobinFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new fixedFluxTransportRobinFvPatchScalarField(*this, iF)
            );
        }


        // Mapping functions

             //- Map (and resize as needed) from self given a mapping object
             virtual void autoMap
             (
                 const fvPatchFieldMapper&
             );

             //- Reverse map the given fvPatchField onto this fvPatchField
             virtual void rmap
             (
                 const fvPatchField<scalar>&,
                 const labelList&
             );


        // Member functions

            // virtual scalarField& RobinK()
            // {
            //     return RobinKeff_;
            // }
            //
            // virtual const scalarField& RobinK() const
            // {
            //     return RobinKeff_;
            // }

            // //- Evaluate the patch field
            // virtual void evaluate
            // (
            //     const Pstream::commsTypes commsType =
            //         Pstream::commsTypes::blocking
            // );

           //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<Field<scalar>> valueInternalCoeffs
            (
                const tmp<scalarField>&
            ) const;

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<Field<scalar>> valueBoundaryCoeffs
            (
                const tmp<scalarField>&
            ) const;

            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<Field<scalar>> gradientInternalCoeffs() const;

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<Field<scalar>> gradientBoundaryCoeffs() const;

            virtual void updateCoeffs();

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
