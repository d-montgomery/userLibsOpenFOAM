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


\*---------------------------------------------------------------------------*/

#include "fixedFluxTransportRobinFvPatchScalarField.H"
// #include "fvPatchFieldMapper.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fixedFluxTransportRobinFvPatchScalarField::
fixedFluxTransportRobinFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    ghostRobinFvPatchScalarField(p, iF),
    phiName_("phi"),
    useRho_(false)
{}


Foam::fixedFluxTransportRobinFvPatchScalarField::
fixedFluxTransportRobinFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    ghostRobinFvPatchScalarField(p, iF, dict),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    useRho_(dict.lookupOrDefault<bool>("useRho", false))
{}


Foam::fixedFluxTransportRobinFvPatchScalarField::
fixedFluxTransportRobinFvPatchScalarField
(
    const fixedFluxTransportRobinFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    ghostRobinFvPatchScalarField(p, iF),
    phiName_(ptf.phiName_),
    useRho_(ptf.useRho_)
{}


Foam::fixedFluxTransportRobinFvPatchScalarField::
fixedFluxTransportRobinFvPatchScalarField
(
    const fixedFluxTransportRobinFvPatchScalarField& ptf
)
:
    ghostRobinFvPatchScalarField(ptf),
    phiName_(ptf.phiName_),
    useRho_(ptf.useRho_)
{}


Foam::fixedFluxTransportRobinFvPatchScalarField::
fixedFluxTransportRobinFvPatchScalarField
(
    const fixedFluxTransportRobinFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    ghostRobinFvPatchScalarField(ptf, iF),
    phiName_(ptf.phiName_),
    useRho_(ptf.useRho_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::fixedFluxTransportRobinFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper&   m
)
{
    ghostRobinFvPatchScalarField::autoMap(m);
//    RobinKeff_.autoMap(m);
}

void Foam::fixedFluxTransportRobinFvPatchScalarField::rmap
(
    const fvPatchField<scalar>& ptf,
    const labelList& addr
)
{
    ghostRobinFvPatchScalarField::rmap(ptf,addr);
    //
    // const fixedFluxTransportRobinFvPatchScalarField& mptf =
    //     refCast<const fixedFluxTransportRobinFvPatchScalarField>(ptf);
    //
    // RobinKeff_.rmap(mptf.RobinKeff_,addr);
}

void Foam::fixedFluxTransportRobinFvPatchScalarField::write(Ostream& os) const
{
    ghostRobinFvPatchScalarField::write(os);
    os.writeEntryIfDifferent<word>("phi", "phi", phiName_);
    os.writeEntry("useRho",useRho_);
    // RobinKeff_.writeEntry("RobinKeff_",os);
}

// void Foam::fixedFluxTransportRobinFvPatchScalarField::evaluate
// (
//     const Pstream::commsTypes commsType
// )
// {
//
//     const scalarField& RobinK = RobinFvPatchScalarField::RobinK();
//
//     const fvsPatchField<scalar>& phip =
//         patch().lookupPatchField<surfaceScalarField, scalar>(phiName_);
//
//     //- Calculate effective Robin coefficient
//     RobinKeff_ = phip/patch().magSf()
//                   + RobinK;
//
//     //- Evaluate Robin boundary condition
//     RobinFvPatchScalarField::evaluate();
//
// }


void Foam::fixedFluxTransportRobinFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
      return;
    }

    const fvsPatchField<scalar>& phip =
        patch().lookupPatchField<surfaceScalarField, scalar>(phiName_);
    const fvPatchField<scalar>& Dp =
        patch().lookupPatchField<volScalarField, scalar>("thermo:D");

     RobinD() = Dp;

    if(useRho_)
    {
        const fvPatchField<scalar>& rhop =
            patch().lookupPatchField<volScalarField, scalar>("rho");

        RobinD() *= rhop;
    }

    RobinK() = phip/patch().magSf();
    RobinF() = 0.;

    //- Evaluate Robin boundary condition
    ghostRobinFvPatchScalarField::updateCoeffs();

}

Foam::tmp<Foam::Field<Foam::scalar>>
Foam::fixedFluxTransportRobinFvPatchScalarField::valueInternalCoeffs
(
    const tmp<scalarField>& f
) const
{
    // if(max(RobinF()) < SMALL)
    // {
    //     return RobinK()*0.;
    // }

    return ghostRobinFvPatchField<scalar>::valueInternalCoeffs(f);
}

Foam::tmp<Foam::Field<Foam::scalar>>
Foam::fixedFluxTransportRobinFvPatchScalarField::valueBoundaryCoeffs
(
    const tmp<scalarField>& f
) const
{
    // if(max(RobinF()) < SMALL)
    // {
    //     return RobinK()*0.;
    // }

    return ghostRobinFvPatchField<scalar>::valueBoundaryCoeffs(f);
}

Foam::tmp<Foam::Field<Foam::scalar>>
Foam::fixedFluxTransportRobinFvPatchScalarField::gradientInternalCoeffs() 
const
{
    // if(max(RobinF()) < SMALL)
    // {
    //     return RobinK()*0.;
    // }

    return ghostRobinFvPatchField<scalar>::gradientInternalCoeffs();
}

Foam::tmp<Foam::Field<Foam::scalar>>
Foam::fixedFluxTransportRobinFvPatchScalarField::gradientBoundaryCoeffs()
const
{
    // if(max(RobinF()) < SMALL)
    // {
    //     return RobinK()*0.;
    // }

    return ghostRobinFvPatchField<scalar>::gradientBoundaryCoeffs();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        fixedFluxTransportRobinFvPatchScalarField
    );
}

// ************************************************************************* //
