/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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
    Foam::ghostRobinFvPatchField

Group
    grpGenericBoundaryConditions

Description
    A Robin boundary condition implemented using ghost nodes for high order accuracy.
    Derivatives at the boundary are represented using central differencing and
    the ghost node is taken at a distance delta from the boundary.

            |
            |
    *       *       *
    C       |P       G (ghost)
            |
            |

    The field at the boundary:

        psi_p = (psi_c + psi_g)/2

    And the Robin condition is

        k(psi_c + psi_g)/2 -  D (psi_g - psi_c)/2dn + f = 0

    So that:

        psi_g =   psi_c (D/dn + k)/(D/dn- k) + 2f/(D/dn-k)

    The normal derivative is:

        dpsi/dn = (psi_g - psi_c)/2dn
                =   psi_c ((D/dn + k)/(D/dn- k) - 1)/2dn + f/(D - dn k)

    And the value is:

        psi_g = psi_c (1 + (D/dn + k)/(D/dn- k))/2 + f/(D/dn + k)


    All boundary conditions can be implemented as derived types of this
    special Robin condition.

    Ghost nodes are virtual and no equation is actually solved on them.
    psi_g is expressed as a function of psi_c and used to compute the boundary value
    and gradient implicitly.

SourceFiles
    ghostRobinFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef ghostRobinFvPatchField_H
#define ghostRobinFvPatchField_H

#include "fvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class ghostRobinFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class ghostRobinFvPatchField
:
    public fvPatchField<Type>
{
    // Private data

    //- Diffusion field
    scalarField RobinD_;

    //- Reaction field
    scalarField RobinK_;

    //- Source term
    Field<Type> RobinF_;


public:

    //- Runtime type information
    TypeName("ghostRobin");


    // Constructors

        //- Construct from patch and internal field
        ghostRobinFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        ghostRobinFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping the given ghostRobinFvPatchField
        //  onto a new patch
        ghostRobinFvPatchField
        (
            const ghostRobinFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        ghostRobinFvPatchField
        (
            const ghostRobinFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type>> clone() const
        {
            return tmp<fvPatchField<Type>>
            (
                new ghostRobinFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        ghostRobinFvPatchField
        (
            const ghostRobinFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type>> clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type>>
            (
                new ghostRobinFvPatchField<Type>(*this, iF)
            );
        }


    // Member functions

    // Return defining fields

        virtual scalarField& RobinD()
        {
            return RobinD_;
        }

        virtual const scalarField& RobinD() const
        {
            return RobinD_;
        }

        virtual scalarField& RobinK()
        {
            return RobinK_;
        }

        virtual const scalarField& RobinK() const
        {
            return RobinK_;
        }

        virtual Field<Type>& RobinF()
        {
            return RobinF_;
        }

        virtual const Field<Type>& RobinF() const
        {
            return RobinF_;
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
                const fvPatchField<Type>&,
                const labelList&
            );


        // Evaluation functions

            //- Return gradient at boundary
            virtual tmp<Field<Type>> snGrad() const;


            //- Evaluate the patch field
            virtual void evaluate
            (
                const Pstream::commsTypes commsType =
                    Pstream::commsTypes::blocking
            );

            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<Field<Type>> valueInternalCoeffs
            (
                const tmp<scalarField>&
            ) const;

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<Field<Type>> valueBoundaryCoeffs
            (
                const tmp<scalarField>&
            ) const;

            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<Field<Type>> gradientInternalCoeffs() const;

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<Field<Type>> gradientBoundaryCoeffs() const;

            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ghostRobinFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
