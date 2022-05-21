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

Class
    Foam::convectiveRadiativeHeatFluxFvPatchScalarField

Description
    Ranbedingung welche einen konvektiven Wärmeübergang mit einem Modell für 
    Wärmestrahlung kombiniert. Für den spezifischen Fall von 2 koaxialen Rohren

    The ambient temperature Ta is specified as a Foam::Function1 of time but
    uniform is space.

Usage
    \table
    Property     | Description          			| Required  | Default value
    mode         | 'convecticeRadiative'			|    yes    |
    Q            | Power [W]          				|     no    |
    q            | Heat flux [W/m^2]           			|     no    |
    h            | Heat transfer coefficient [W/m^2/K] 		|    yes    |
    Ta           | Ambient temperature [K]     			|    yes    |
    relaxation   | Relaxation for the wall temperature 		|     no    |      1
    emissivityNozzle   | Surface emissivity Nozzle		|    yes    |
    emissivityFilament | Surface emissivity Filament 		|    yes    |
    kappaMethod  | Inherited from temperatureCoupledBase 	| inherited |
    kappa        | Inherited from temperatureCoupledBase 	| inherited |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            	convectiveRadiativeHeatFlux;
        mode           		convectiveRadiative;
        Ta              	constant 510.0;
        h               	uniform 10.0;
	emissivityNozzle	0.9;
	emissivityFilament	0.6;
	dNozzle			3.2;
	dFilament		1.75;
        kappaMethod     	fluidThermo;
        value           	$internalField;		//initial Temperature value
    }
    \endverbatim

See also
    Foam::temperatureCoupledBase
    Foam::mixedFvPatchScalarField

SourceFiles
    convectiveRadiativeHeatFluxFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef convectiveRadiativeHeatFluxFvPatchScalarField_H
#define convectiveRadiativeHeatFluxFvPatchScalarField_H

#include "mixedFvPatchFields.H"
#include "temperatureCoupledBase.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
      Class convectiveRadiativeHeatFluxFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class convectiveRadiativeHeatFluxFvPatchScalarField
:
    public mixedFvPatchScalarField,
    public temperatureCoupledBase
{
public:

    // Public data

        //- Operation mode enumeration
        enum operationMode
        {
            convectiveRadiative
        };

        static const NamedEnum<operationMode, 1> operationModeNames;		// 1 Mode


private:

    // Private Data

        //- Operation mode
        operationMode mode_;

        //- Heat power [W]
        scalar Q_;

        //- Heat flux [W/m^2]
        scalarField q_;

        //- Heat transfer coefficient [W/m^2K]
        scalarField h_;

        //- Ambient temperature [K]
        autoPtr<Function1<scalar>> Ta_;

        //- Optional surface emissivity for radiative transfer to ambient
        scalar emissivityNozzle_;

        //- Optional surface emissivity for radiative transfer to ambient
        scalar emissivityFilament_;

        //- Optional surface emissivity for radiative transfer to ambient
        scalar dNozzle_;

        //- Optional surface emissivity for radiative transfer to ambient
        scalar dFilament_;

        //- Relaxation for the wall temperature (thermal inertia)
        scalar relaxation_;


public:

    //- Runtime type information
    TypeName("convectiveRadiativeHeatFlux");


    // Constructors

        //- Construct from patch and internal field
        convectiveRadiativeHeatFluxFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        convectiveRadiativeHeatFluxFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  convectiveRadiativeHeatFluxFvPatchScalarField
        //  onto a new patch
        convectiveRadiativeHeatFluxFvPatchScalarField
        (
            const convectiveRadiativeHeatFluxFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Copy constructor
        convectiveRadiativeHeatFluxFvPatchScalarField
        (
            const convectiveRadiativeHeatFluxFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new convectiveRadiativeHeatFluxFvPatchScalarField(*this)
            );
        }

        //- Copy constructor setting internal field reference
        convectiveRadiativeHeatFluxFvPatchScalarField
        (
            const convectiveRadiativeHeatFluxFvPatchScalarField&,
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
                new convectiveRadiativeHeatFluxFvPatchScalarField(*this, iF)
            );
        }


    // Member Functions

        // Access

            //- Allow manipulation of the boundary values
            virtual bool fixesValue() const
            {
                return false;
            }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            //  Used to update fields following mesh topology change
            virtual void autoMap(const fvPatchFieldMapper&);

            //- Reverse map the given fvPatchField onto this fvPatchField
            //  Used to reconstruct fields
            virtual void rmap(const fvPatchScalarField&, const labelList&);


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        // I-O

            //- Write
            void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
