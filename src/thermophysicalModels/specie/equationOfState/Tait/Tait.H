/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2019 OpenFOAM Foundation
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
    Foam::Tait

Description
    Vereinfachte Tait Zustandsgleichung (EOS) für amorphe Polymere.
    Vereinfachungen:
	Nur für amorphe Polymere, Modell könnte durch einbinden der b7-, b8- und b9- Parameter auf Teilkristalline Kunststoffe erweitert werden 
	Nur Temperaturabhängig, da die Druckabhängigkeit in der Arbeit keinen großen Einfluss hatte
	Druckabhängigkeit allerdings implementiert! Es müssen nur die entsprechend gekennzeichneten auskommentierten Stellen im Code wieder eingebunden werden

SourceFiles
    TaitI.H
    Tait.C

\*---------------------------------------------------------------------------*/

#ifndef Tait_H
#define Tait_H

#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Specie> class Tait;

template<class Specie>
inline Tait<Specie> operator+
(
    const Tait<Specie>&,
    const Tait<Specie>&
);

template<class Specie>
inline Tait<Specie> operator*
(
    const scalar,
    const Tait<Specie>&
);

template<class Specie>
inline Tait<Specie> operator==
(
    const Tait<Specie>&,
    const Tait<Specie>&
);

template<class Specie>
Ostream& operator<<
(
    Ostream&,
    const Tait<Specie>&
);



/*---------------------------------------------------------------------------*\
                           Class Tait Declaration
\*---------------------------------------------------------------------------*/

template<class Specie>
class Tait
:
    public Specie
{
    // Private Data

        // Tait Koeffizient b5
        scalar b5_;

        // Tait Koeffizient b6
        scalar b6_;

        // Tait Koeffizient b1m
        scalar b1m_;

        // Tait Koeffizient b2m
        scalar b2m_;

        // Tait Koeffizient b3m
        scalar b3m_;

        // Tait Koeffizient b4m
        scalar b4m_;

        // Tait Koeffizient b1s
        scalar b1s_;

        // Tait Koeffizient b2s
        scalar b2s_;

        // Tait Koeffizient b3s
        scalar b3s_;

        // Tait Koeffizient b4s
        scalar b4s_;

public:

    // Constructors

        //- Construct from components
        inline Tait
        (
            const Specie& sp,
            const scalar& b5,
            const scalar& b6,
            const scalar& b1m,
            const scalar& b2m,
            const scalar& b3m,
            const scalar& b4m,
            const scalar& b1s,
            const scalar& b2s,
            const scalar& b3s,
            const scalar& b4s
        );

        //- Construct from dictionary
        Tait(const dictionary& dict);

        //- Construct as named copy
        inline Tait(const word& name, const Tait&);

        //- Construct and return a clone
        inline autoPtr<Tait> clone() const;

        // Selector from dictionary
        inline static autoPtr<Tait> New
        (
            const dictionary& dict
        );


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "Tait<" + word(Specie::typeName_()) + '>';
        }

        // Fundamental properties


            //- Is the equation of state is incompressible i.e. rho != f(p)
            static const bool incompressible = true;

            //- Is the equation of state is isochoric i.e. rho = const
            static const bool isochoric = false;

            //- Return density [kg/m^3]
            inline scalar rho(scalar p, scalar T) const;

            //- Return enthalpy departure [J/kg]
            inline scalar H(const scalar p, const scalar T) const;

            //- Return Cp departure [J/(kg K]
            inline scalar Cp(scalar p, scalar T) const;

            //- Return internal energy departure [J/kg]
            inline scalar E(const scalar p, const scalar T) const;

            //- Return Cv departure [J/(kg K]
            inline scalar Cv(scalar p, scalar T) const;

            //- Return entropy [J/kg/K]
            inline scalar S(const scalar p, const scalar T) const;

            //- Return compressibility rho/p [s^2/m^2]
            inline scalar psi(scalar p, scalar T) const;

            //- Return compression factor []
            inline scalar Z(scalar p, scalar T) const;

            //- Return v_0(T) for Tait-equation
            inline scalar v_0(scalar p, scalar T) const;

            //- Return B(T) for Tait-equation  			//____________einbinden für Druckabhängigkeit
            //inline scalar B(scalar p, scalar T) const;

            //- Return (Cp - Cv) [J/(kg K]
            inline scalar CpMCv(scalar p, scalar T) const;


    // IO

            //- Write to Ostream
            void write(Ostream& os) const;

    // Member Operators

        inline void operator+=(const Tait&);
        inline void operator*=(const scalar);


    // Friend operators

        friend Tait operator+ <Specie>
        (
            const Tait&,
            const Tait&
        );

        friend Tait operator* <Specie>
        (
            const scalar s,
            const Tait&
        );

        friend Tait operator== <Specie>
        (
            const Tait&,
            const Tait&
        );


    // Ostream Operator

        friend Ostream& operator<< <Specie>
        (
            Ostream&,
            const Tait&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "TaitI.H"

#ifdef NoRepository
    #include "Tait.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
