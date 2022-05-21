/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2019 OpenFOAM Foundation
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
    Foam::KarrenbergWLFTransport

Description

    Modifiziertes WLF-Modell basierend auf dem Ansatz nach Gregor Karrenberg mit eigener Modellierung
    Siehe schriftliche Ausarbeitung Masterarbeit Tim Ebinger für die entsprechende Gleichung

    References:
    \verbatim
        Williams, M. L., Landel, R. F., & Ferry, J. D. (1955).
        The temperature dependence of relaxation mechanisms
        in amorphous polymers and other glass-forming liquids.
        Journal of the American Chemical society, 77(14), 3701-3707.
    \endverbatim

SourceFiles
    KarrenbergWLFTransportI.H
    KarrenbergWLFTransport.C

\*---------------------------------------------------------------------------*/

#ifndef KarrenbergWLFTransport_H
#define KarrenbergWLFTransport_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Thermo> class KarrenbergWLFTransport;

template<class Thermo>
inline KarrenbergWLFTransport<Thermo> operator+
(
    const KarrenbergWLFTransport<Thermo>&,
    const KarrenbergWLFTransport<Thermo>&
);

template<class Thermo>
inline KarrenbergWLFTransport<Thermo> operator*
(
    const scalar,
    const KarrenbergWLFTransport<Thermo>&
);

template<class Thermo>
Ostream& operator<<
(
    Ostream&,
    const KarrenbergWLFTransport<Thermo>&
);


/*---------------------------------------------------------------------------*\
                    Class KarrenbergWLFTransport Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class KarrenbergWLFTransport
:
    public Thermo
{
    // Private Data

        //- Dynamic viscosity at the reference temperature [Pa.s]
        scalar mu0_;

        //- Reference temperature [T]
        scalar Tr_;

        //-  Melt-Fraction Koeffizient für Karrenberg-Modellierung Tm
        scalar Tm_;

        //-  Melt-Fraction Koeffizient für Karrenberg-Modellierung TS
        scalar TS_;

        //- Melt-Fraction Koeffizient für Karrenberg-Modellierung MS
        scalar MS_;

        //- WLF coefficient 1 []
        scalar C1_;

        //- WLF coefficient 2 [T]
        scalar C2_;



    // Private Member Functions

        //- Read coefficient from dictionary
        scalar readCoeff(const word& coeffName, const dictionary& dict);


public:

    // Constructors

        //- Construct as named copy
        inline KarrenbergWLFTransport(const word&, const KarrenbergWLFTransport&);

        //- Construct from dictionary
        KarrenbergWLFTransport(const dictionary& dict);

        //- Construct and return a clone
        inline autoPtr<KarrenbergWLFTransport> clone() const;

        // Selector from dictionary
        inline static autoPtr<KarrenbergWLFTransport> New(const dictionary& dict);


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "KarrenbergWLF<" + Thermo::typeName() + '>';
        }

        //- Dynamic viscosity [kg/m/s]
        inline scalar mu(const scalar p, const scalar T) const;

        //- Melt-Fraction
        inline scalar MF(const scalar p, const scalar T) const;

        //- Thermal conductivity [W/m/K]
        inline scalar kappa(const scalar p, const scalar T) const;

        //- Thermal diffusivity of enthalpy [kg/m/s]
        inline scalar alphah(const scalar p, const scalar T) const;

        // Species diffusivity
        // inline scalar D(const scalar p, const scalar T) const;

        //- Write to Ostream
        void write(Ostream& os) const;


    // Member Operators

        inline void operator+=(const KarrenbergWLFTransport&);
        inline void operator*=(const scalar);


    // Friend operators

        friend KarrenbergWLFTransport operator+ <Thermo>
        (
            const KarrenbergWLFTransport&,
            const KarrenbergWLFTransport&
        );

        friend KarrenbergWLFTransport operator* <Thermo>
        (
            const scalar,
            const KarrenbergWLFTransport&
        );


    // Ostream Operator

        friend Ostream& operator<< <Thermo>
        (
            Ostream&,
            const KarrenbergWLFTransport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "KarrenbergWLFTransportI.H"

#ifdef NoRepository
    #include "KarrenbergWLFTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
