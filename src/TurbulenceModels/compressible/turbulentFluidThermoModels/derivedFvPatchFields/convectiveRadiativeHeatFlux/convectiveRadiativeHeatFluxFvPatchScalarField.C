/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

#include "convectiveRadiativeHeatFluxFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "physicoChemicalConstants.H"

using Foam::constant::physicoChemical::sigma;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam	//_________________Operation Modes definiert als Enum --> Abfrage enum in Switch-Case Anweisungen
{
    template<>
    const char*
    NamedEnum
    <
        convectiveRadiativeHeatFluxFvPatchScalarField::operationMode,
       1										// number of Modes
    >::names[] =
    {
        "convectiveRadiative"								// Bezeichnung Modes
    };
}

const Foam::NamedEnum
<
    Foam::convectiveRadiativeHeatFluxFvPatchScalarField::operationMode,
   1											// number of Modes
> Foam::convectiveRadiativeHeatFluxFvPatchScalarField::operationModeNames;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//___________________________________________________________________________________Initialisierung 
Foam::convectiveRadiativeHeatFluxFvPatchScalarField::
convectiveRadiativeHeatFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    temperatureCoupledBase(patch(), dict),
    mode_(operationModeNames.read(dict.lookup("mode"))),
    Q_(0),
    Ta_(),
    emissivityNozzle_(),		
    emissivityFilament_(),
    dNozzle_(),
    dFilament_(),		
    relaxation_(dict.lookupOrDefault<scalar>("relaxation", 1))
{
    switch (mode_)
    {
        case convectiveRadiative:
        {
            h_ = scalarField("h", dict, p.size());
            Ta_ = Function1<scalar>::New("Ta", dict);
	    
	    if (dict.found("emissivityNozzle"))
            {
                dict.lookup("emissivityNozzle") >> emissivityNozzle_;
                dict.lookup("emissivityFilament") >> emissivityFilament_;
                dict.lookup("dNozzle") >> dNozzle_;
                dict.lookup("dFilament") >> dFilament_;
            }
	    else
            {
                Info<< "Bitte Parameter X angeben" << endl;
            }
            break;
        }
    }

    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));

    if (dict.found("refValue"))
    {
        // Full restart
        refValue() = scalarField("refValue", dict, p.size());
        refGrad() = scalarField("refGradient", dict, p.size());
        valueFraction() = scalarField("valueFraction", dict, p.size());
    }
    else
    {
        // Start from user entered data. Assume fixedValue.
        refValue() = *this;
        refGrad() = 0;
        valueFraction() = 1;
    }
}
//___________________________________________________________________________________bis hier Initialisierung

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//___________________________________________________________________________________mappen der Wärmeübergangskoeffizienten auf Feld
void Foam::convectiveRadiativeHeatFluxFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mixedFvPatchScalarField::autoMap(m);

    switch (mode_)
    {
        case convectiveRadiative:
        {
            m(h_, h_);

            break;
        }
    }
}

void Foam::convectiveRadiativeHeatFluxFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    mixedFvPatchScalarField::rmap(ptf, addr);

    const convectiveRadiativeHeatFluxFvPatchScalarField& tiptf =
        refCast<const convectiveRadiativeHeatFluxFvPatchScalarField>(ptf);

    switch (mode_)
    {
        case convectiveRadiative:
        {
            h_.rmap(tiptf.h_, addr);

            break;
        }
    }
}
//___________________________________________________________________________________bis hier mappen
//___________________________________________________________________________________ab hier Berechnung
void Foam::convectiveRadiativeHeatFluxFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const scalarField& Tp(*this);

    // Store current valueFraction and refValue for relaxation
    const scalarField valueFraction0(valueFraction());
    const scalarField refValue0(refValue());

    switch (mode_)
    {
        case convectiveRadiative:
        {
            scalarField hp(h_);

            const scalar Ta = Ta_->value(this->db().time().timeOutputValue());
            scalarField hpTa(hp*Ta);

            const scalarField kappaDeltaCoeffs
            (
                this->kappa(Tp)*patch().deltaCoeffs()
            );

	    const scalar A1_ = dFilament_;		//Übergabe des Filamentdurchmessers
	    const scalar A2_ = dNozzle_;		//Übergabe des Düsendurchmessers
	    const scalar C12_ = sigma.value()/((1/emissivityFilament_)+(A1_/A2_)*((1/emissivityNozzle_)-1));	//Berechnung des Strahlungsaustauschkoeffizienten
	    refGrad() = (C12_*(pow4(Ta)-pow4(Tp)))/kappa(Tp);		// Berechnung des Strahlungsanteils	(durch Kappa teilen da gemäß Definition gradT = q/Kappa)

            forAll(Tp, i)						//für mixedFvPatchField.H bzw. .C in src/finiteVolume/fields/fvPatchFields/basic/mixed
            {
                    refValue()[i] = (hpTa[i])/hp[i];
                    valueFraction()[i] = hp[i]/(hp[i] + kappaDeltaCoeffs[i]);
            }

            break;
        }
    }
//Relaxation
    valueFraction() =
        relaxation_*valueFraction()
      + (1 - relaxation_)*valueFraction0;
    refValue() = relaxation_*refValue() + (1 - relaxation_)*refValue0;
//Ende Relaxation

    mixedFvPatchScalarField::updateCoeffs();

    if (debug)		//für debugging
    {
        const scalar Q = gSum(kappa(Tp)*patch().magSf()*snGrad());
	const scalar A1_ = dFilament_;
	const scalar A2_ = dNozzle_;
	const scalar C12_ = sigma.value()/((1/emissivityFilament_)+(A1_/A2_)*((1/emissivityNozzle_)-1));

        Info<< patch().boundaryMesh().mesh().name() << ':'	//region0
            << patch().name() << ':'				//wall_Top
            << this->internalField().name() << " :"		//heat Transfer Rate
            << " heat transfer rate:" << Q		//
            << "  weitere Parameter"
            << " C12" << C12_				//Strahlungskonstante
            << " Tp" << Tp				//Temperatur Wand
            << " min:" << gMin(*this)			//aktuelles Tmin
            << " max:" << gMax(*this)			//aktuelles Tmax
            << " kappa(Tp)" << kappa(Tp)		//
            << " patch().magSf()" << patch().magSf()	//	=   Zellfläche Wand-Patches
            << " snGrad:" << snGrad()			//surfaceNormalGRadient	
            << endl;
    }
}
//___________________________________________________________________________________bis hier Berechnung
//___________________________________________________________________________________ab hier Ausgabe
void Foam::convectiveRadiativeHeatFluxFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);

    writeEntry(os, "mode", operationModeNames[mode_]);
    temperatureCoupledBase::write(os);

    switch (mode_)
    {
        case convectiveRadiative:
        {
            writeEntry(os, "h", h_);
            writeEntry(os, Ta_());
	    writeEntry(os, "emissivityNozzle", emissivityNozzle_);
	    writeEntry(os, "emissivityFilament", emissivityFilament_);
	    writeEntry(os, "dNozzle", dNozzle_);
	    writeEntry(os, "dFilament", dFilament_);
            break;
        }
    }

    writeEntry(os, "refValue", refValue());
    writeEntry(os, "refGradient", refGrad());
    writeEntry(os, "valueFraction", valueFraction());
    writeEntry(os, "value", *this);
}
//___________________________________________________________________________________bis hier Ausgabe

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        convectiveRadiativeHeatFluxFvPatchScalarField
    );
}

// ************************************************************************* //
