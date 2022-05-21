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
    Foam::hTaitThermo

Description
    Berücksichtigung einer temperaturabhängigen Enthalpie
    Berechnet aus Cp(T) mittels MatLAB --> Integration über Temperaturbereich --> Ausgabe eines Arrays (siehe unten)

SourceFiles
    hTaitThermoI.H
    hTaitThermo.C

\*---------------------------------------------------------------------------*/

#ifndef hTaitThermo_H
#define hTaitThermo_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class EquationOfState> class hTaitThermo;

template<class EquationOfState>
inline hTaitThermo<EquationOfState> operator+
(
    const hTaitThermo<EquationOfState>&,
    const hTaitThermo<EquationOfState>&
);

template<class EquationOfState>
inline hTaitThermo<EquationOfState> operator*
(
    const scalar,
    const hTaitThermo<EquationOfState>&
);

template<class EquationOfState>
inline hTaitThermo<EquationOfState> operator==
(
    const hTaitThermo<EquationOfState>&,
    const hTaitThermo<EquationOfState>&
);

template<class EquationOfState>
Ostream& operator<<
(
    Ostream&,
    const hTaitThermo<EquationOfState>&
);


/*---------------------------------------------------------------------------*\
                           Class hTaitThermo Declaration
\*---------------------------------------------------------------------------*/

template<class EquationOfState>
class hTaitThermo
:
    public EquationOfState
{
    // Private Data

        scalar Cp_;


    // Private Member Functions

        //- Construct from components
        inline hTaitThermo
        (
            const EquationOfState& st,
            const scalar cp
        );


public:

    // Constructors

        //- Construct from dictionary
        hTaitThermo(const dictionary& dict);

        //- Construct as named copy
        inline hTaitThermo(const word&, const hTaitThermo&);

        //- Construct and return a clone
        inline autoPtr<hTaitThermo> clone() const;

        //- Selector from dictionary
        inline static autoPtr<hTaitThermo> New(const dictionary& dict);


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "hTait<" + EquationOfState::typeName() + '>';
        }

        //- Limit the temperature to be in the range Tlow_ to Thigh_
        inline scalar limit(const scalar T) const;


        // Fundamental properties

            //- Heat capacity at constant pressure [J/kg/K]
            inline scalar Cp(const scalar p, const scalar T) const;

            //- Absolute Enthalpy [J/kg]
            inline scalar Ha(const scalar p, const scalar T) const;

            //- Sensible enthalpy [J/kg]
            inline scalar Hs(const scalar p, const scalar T) const;

            //- Chemical enthalpy [J/kg]
            inline scalar Hc() const;

            //- Entropy [J/kg/K]
            inline scalar S(const scalar p, const scalar T) const;

            #include "HtoEthermo.H"


        // Derivative term used for Jacobian

            //- Derivative of Gibbs free energy w.r.t. temperature
            inline scalar dGdT(const scalar p, const scalar T) const;

            //- Temperature derivative of heat capacity at constant pressure
            inline scalar dCpdT(const scalar p, const scalar T) const;


        // I-O

            //- Write to Ostream
            void write(Ostream& os) const;


    // Member Operators

        inline void operator+=(const hTaitThermo&);


    // Friend operators

        friend hTaitThermo operator+ <EquationOfState>
        (
            const hTaitThermo&,
            const hTaitThermo&
        );

        friend hTaitThermo operator* <EquationOfState>
        (
            const scalar,
            const hTaitThermo&
        );

        friend hTaitThermo operator== <EquationOfState>
        (
            const hTaitThermo&,
            const hTaitThermo&
        );


    // IOstream Operators

        friend Ostream& operator<< <EquationOfState>
        (
            Ostream&,
            const hTaitThermo&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "hTaitThermoI.H"

#ifdef NoRepository
    #include "hTaitThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
