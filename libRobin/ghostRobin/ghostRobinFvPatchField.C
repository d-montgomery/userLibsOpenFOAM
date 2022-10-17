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

\*---------------------------------------------------------------------------*/

#include "ghostRobinFvPatchField.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::ghostRobinFvPatchField<Type>::ghostRobinFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(p, iF),
    RobinD_(p.size()),
    RobinK_(p.size()),
    RobinF_(p.size())
{}


template<class Type>
Foam::ghostRobinFvPatchField<Type>::ghostRobinFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fvPatchField<Type>(p, iF, dict, false),
    RobinD_(p.size()),
    RobinK_(p.size()),
    RobinF_(p.size())
{
    if (dict.found("RobinD"))
    {
      RobinD_ = scalarField("RobinD",dict,p.size());
    }
    else
    {
      RobinD_ = scalarField(p.size(),1.0);
    }
    if (dict.found("RobinK"))
    {
      RobinK_ = scalarField("RobinK",dict,p.size());
    }
    else
    {
      RobinK_ = scalarField(p.size(),0.0);
    }
    if (dict.found("ghostRobinF"))
    {
      RobinF_ = Field<Type>("ghostRobinF",dict,p.size());
    }
    else
    {
      RobinF_ = Field<Type>(p.size(),pTraits<Type>::zero);
    }

    evaluate();
}


template<class Type>
Foam::ghostRobinFvPatchField<Type>::ghostRobinFvPatchField
(
    const ghostRobinFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fvPatchField<Type>(ptf, p, iF, mapper),
    RobinD_(ptf.RobinD_, mapper),
    RobinK_(ptf.RobinK_, mapper),
    RobinF_(ptf.RobinF_, mapper)
{
    if (notNull(iF) && mapper.hasUnmapped())
    {
        WarningInFunction
            << "On field " << iF.name() << " patch " << p.name()
            << " patchField " << this->type()
            << " : mapper does not map all values." << nl
            << "    To avoid this warning fully specify the mapping in derived"
            << " patch fields." << endl;
    }
}


template<class Type>
Foam::ghostRobinFvPatchField<Type>::ghostRobinFvPatchField
(
    const ghostRobinFvPatchField<Type>& ptf
)
:
    fvPatchField<Type>(ptf),
    RobinD_(ptf.RobinD_),
    RobinK_(ptf.RobinK_),
    RobinF_(ptf.RobinF_)
{}


template<class Type>
Foam::ghostRobinFvPatchField<Type>::ghostRobinFvPatchField
(
    const ghostRobinFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(ptf, iF),
    RobinD_(ptf.RobinD_),
    RobinK_(ptf.RobinK_),
    RobinF_(ptf.RobinF_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::ghostRobinFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fvPatchField<Type>::autoMap(m);
    RobinD_.autoMap(m);
    RobinK_.autoMap(m);
    RobinF_.autoMap(m);
}


template<class Type>
void Foam::ghostRobinFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    fvPatchField<Type>::rmap(ptf, addr);

    const ghostRobinFvPatchField<Type>& fgptf =
        refCast<const ghostRobinFvPatchField<Type>>(ptf);

        RobinD_.rmap(fgptf .RobinD_, addr);
        RobinK_.rmap(fgptf .RobinK_, addr);
        RobinF_.rmap(fgptf .RobinF_, addr);
}


template<class Type>
void Foam::ghostRobinFvPatchField<Type>::evaluate(const Pstream::commsTypes)
{

    if (!this->updated())
    {
        this->updateCoeffs();
    }

    scalarField coeffP(RobinD() * this->patch().deltaCoeffs() + RobinK());
    scalarField coeffN(RobinD() * this->patch().deltaCoeffs() - RobinK());



    Field<Type>::operator=
    (
        // (
        //   this->patchInternalField() * RobinD() * this->patch().deltaCoeffs()
        //   + RobinF()
        // )
        // /
        // (
        //    RobinD() * this->patch().deltaCoeffs()
        //   - RobinK()
        // )
        this->patchInternalField()*(1. + coeffP/coeffN)*0.5
        + RobinF()/coeffN
    );

    fvPatchField<Type>::evaluate();

    //- Give error on Robin (homogeneous)
    // const scalarField err
    // (
    //     mag(
    //         (this->snGrad()*RobinD_) - (RobinK_*(*this))
    //     )/((RobinK_*mag(*this))+SMALL)
    // );
    
    // Pout<< "Max/average error on Robin : " << max(mag(err)) 
    //     << "    " << sum(mag(err))/err.size() << endl;   

 //   Pout << "Coeff : " << max((1. + coeffP/coeffN))*0.5 << endl;

}

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::ghostRobinFvPatchField<Type>::snGrad() const
{

//    scalarField coeffP(RobinD() * this->patch().deltaCoeffs() + RobinK());
//    scalarField coeffN(RobinD() * this->patch().deltaCoeffs() - RobinK());

    return
    (
        ((*this) - this->patchInternalField())*this->patch().deltaCoeffs()
    );
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::ghostRobinFvPatchField<Type>::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    scalarField coeffP(RobinD() * this->patch().deltaCoeffs() + RobinK());
    scalarField coeffN(RobinD() * this->patch().deltaCoeffs() - RobinK());

    return Type(pTraits<Type>::one) *
    (
        (1. + coeffP/coeffN)*0.5
    );

}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::ghostRobinFvPatchField<Type>::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    scalarField coeffN(RobinD() * this->patch().deltaCoeffs() - RobinK());

    return
    (
        RobinF()/coeffN
    );
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::ghostRobinFvPatchField<Type>::gradientInternalCoeffs() const
{
    scalarField coeffP(RobinD() * this->patch().deltaCoeffs() + RobinK());
    scalarField coeffN(RobinD() * this->patch().deltaCoeffs() - RobinK());

    return Type(pTraits<Type>::one) *
    (
        (coeffP/coeffN - 1.)*(0.5*this->patch().deltaCoeffs())
    );
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::ghostRobinFvPatchField<Type>::gradientBoundaryCoeffs() const
{
    scalarField coeffN(RobinD() * this->patch().deltaCoeffs() - RobinK());
    return
    (
       (RobinF()/coeffN)*(0.5*this->patch().deltaCoeffs())
    );
}


template<class Type>
void Foam::ghostRobinFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    this->writeEntry("value", os);
    RobinD_.writeEntry("RobinD", os);
    RobinK_.writeEntry("RobinK", os);
    RobinF_.writeEntry("ghostRobinF", os);

}

template<class Type>
void Foam::ghostRobinFvPatchField<Type>::updateCoeffs()
 {

     if (this->updated())
     {
         return;
     }

     fvPatchField<Type>::updateCoeffs();
 }

// ************************************************************************* //
