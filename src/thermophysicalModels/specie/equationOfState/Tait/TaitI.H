/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2018 OpenFOAM Foundation
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

#include "Tait.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Specie>
inline Foam::Tait<Specie>::Tait
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
)
:
    Specie(sp),
    b5_(b5),
    b6_(b6),
    b1m_(b1m),
    b2m_(b2m),
    b3m_(b3m),
    b4m_(b4m),
    b1s_(b1s),
    b2s_(b2s),
    b3s_(b3s),
    b4s_(b4s)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
inline Foam::Tait<Specie>::Tait
(
    const word& name,
    const Tait& pg
)
:
    Specie(name, pg),
    b5_(pg.b5_),
    b6_(pg.b6_),
    b1m_(pg.b1m_),
    b2m_(pg.b2m_),
    b3m_(pg.b3m_),
    b4m_(pg.b4m_),
    b1s_(pg.b1s_),
    b2s_(pg.b2s_),
    b3s_(pg.b3s_),
    b4s_(pg.b4s_)
{}


template<class Specie>
inline Foam::autoPtr<Foam::Tait <Specie>>
Foam::Tait<Specie>::clone() const
{
    return autoPtr<Tait<Specie>>
    (
        new Tait<Specie>(*this)
    );
}


template<class Specie>
inline Foam::autoPtr<Foam::Tait<Specie>>
Foam::Tait<Specie>::New
(
    const dictionary& dict
)
{
    return autoPtr<Tait<Specie>>
    (
        new Tait<Specie>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>			//_______Berechnung der Dichte rho mit einem reduzierten Tait-Modell
inline Foam::scalar Foam::Tait<Specie>::rho(scalar p, scalar T) const
{
	const scalar C = 0.0894;

	//const scalar B = this->B(p, T);				//____________einbinden für Druckabhängigkeit
	//return 1/(v_0*(1-C*log(1+(p/B)))); 				//____________einbinden für Druckabhängigkeit, entsprechend untenstehendes auskommentieren

	return 1/((this->v_0(p, T))*(1-C*log(1+(0.0008)))); 		//ohne Druckabhängigkeit, 0.0008 ist Mittelwert und ohne v_T(p,T) da keine kristallinen Polymere

}


template<class Specie>
inline Foam::scalar Foam::Tait<Specie>::H(scalar p, scalar T) const	//_______Berechnung des Beitrags zur Enthalpie, vgl. 1. Hauptsatz der Thermodynamik H = U + p/rho
{
        //Info<< "H" << endl;	//Prüfung ob der Solver diese Stelle im Code aufruft
	return p/(this->rho(p, T));
}

template<class Specie>
inline Foam::scalar Foam::Tait<Specie>::Cp(scalar p, scalar T) const	//_______Lineare Interpolation für untenstehenden Array
{									//_______Messdaten wurden mittels Origin-Spline-Funktion auf Datensatz mit delta_T = 1 K erweitert
     	//Info<< "Cp" << endl;	//Prüfung ob der Solver diese Stelle im Code aufruft
     return 0;		//Cp(T) könnte auch hier implementiert werden statt in hTaitThermoI.H, dafür dort unter Cp : "return EquationOfState::Cp(p, T);"
}


template<class Specie>
inline Foam::scalar Foam::Tait<Specie>::E(scalar p, scalar T) const	//_______Berechnung E________wird bei Verwendung von eConst o.ä. vom Solver aufgerufen
{
        //Info<< "Läuft bis hier 4" << endl;
	return 0;
}


template<class Specie>
inline Foam::scalar Foam::Tait<Specie>::Cv(scalar p, scalar T) const	//_______Berechnung Cv_______wird bei Verwendung von eConst o.ä. vom Solver aufgerufen
{
	return 0;
}


template<class Specie>
inline Foam::scalar Foam::Tait<Specie>::S				//_______Berechnung S_______wird nicht aufgerufen
(
    scalar p,
    scalar T
) const
{
	return 0;
}


template<class Specie>
inline Foam::scalar Foam::Tait<Specie>::psi				//_______Berechnung psi_______wird vom Solver aufgerufen
(
    scalar p,
    scalar T
) const
{
	return 0;
}




template<class Specie>
inline Foam::scalar Foam::Tait<Specie>::Z				//_______Berechnung Z_______wird nicht aufgerufen
(
    scalar p,
    scalar T
) const
{
	return 0;
}

template<class Specie>
inline Foam::scalar Foam::Tait<Specie>::v_0				//_______Berechnung v_0(T) für Tait-Modell
(
    scalar p,
    scalar T
) const
{
    //const scalar Tt = b5_+b6_*p;		//____________einbinden für Druckabhängigkeit, entsprechend untenstehendes auskommentieren
    const scalar Tt = b5_+b6_*1e5;		//konstant 1 bar = 1e5 Pascal
    const scalar Tstrich = T-b5_;

    if (T < Tt)		// solid
    {        					//Info<< "Läuft bis hier 9 " << b1s_+b2s_*Tstrich << endl;		//für Debugging
        return b1s_+b2s_*Tstrich;
    }
    else		// melt
    {       					// Info<< "Läuft bis hier 10 " << b1m_+b2m_*Tstrich << endl;		//für Debugging
        return b1m_+b2m_*Tstrich;
    }
}
/*
template<class Specie>
inline Foam::scalar Foam::Tait<Specie>::B				//_______Berechnung B(T) für Tait-Modell____________einbinden für Druckabhängigkeit
(
    scalar p,
    scalar T
) const
{
    const scalar Tt = b5_+b6_*p;
    const scalar Tstrich = T-b5_;

    if (T < Tt)		// solid
    {        
						//Info<< "Läuft bis hier 11 " << b3s_*exp(-b4s_*Tstrich) << endl;	//für Debugging
        return b3s_*exp(-b4s_*Tstrich);
    }
    else		// melt
    {						//Info<< "Läuft bis hier 12 " << b3m_*exp(-b4m_*Tstrich) << endl;	//für Debugging
        return b3m_*exp(-b4m_*Tstrich);
    }
}
*/
template<class Specie>
inline Foam::scalar Foam::Tait<Specie>::CpMCv				//_______Berechnung CpMCv_______wird vom Solver aufgerufen
(
    scalar p,
    scalar T
) const
{
	return 0;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Specie>
inline void Foam::Tait<Specie>::operator+=
(
    const Tait<Specie>& pg
)
{
    scalar Y1 = this->Y();
    Specie::operator+=(pg);

    if (mag(this->Y()) > small)			//unrelevant da keine Mischung
    {
        Y1 /= this->Y();
        const scalar Y2 = pg.Y()/this->Y();
    }
}


template<class Specie>
inline void Foam::Tait<Specie>::operator*=(const scalar s)
{
     Specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


template<class Specie>
Foam::Tait<Specie> Foam::operator+
(
    const Tait<Specie>& pg1,
    const Tait<Specie>& pg2
)
{
    Specie sp
    (
        static_cast<const Specie&>(pg1)
      + static_cast<const Specie&>(pg2)
    );

    if (mag(sp.Y()) < small)
    {
        return Tait<Specie>
        (
            sp,
            pg1.b5_,
            pg1.b6_,
            pg1.b1m_,
            pg1.b2m_,
            pg1.b3m_,
            pg1.b4m_,
            pg1.b1s_,
            pg1.b2s_,
            pg1.b3s_,
            pg1.b4s_
        );
    }
    else	//unrelevant da keine Mischung
    {
        const scalar Y1 = pg1.Y()/sp.Y();
        const scalar Y2 = pg2.Y()/sp.Y();

        return Tait<Specie>
        (
            sp
        );
    }
}


template<class Specie>
Foam::Tait<Specie> Foam::operator*
(
    const scalar s,
    const Tait<Specie>& pg
)
{
    return Tait<Specie>
    (
        s*static_cast<const Specie&>(pg),
        pg.b5_,
        pg.b6_,
        pg.b1m_,
        pg.b2m_,
        pg.b3m_,
        pg.b4m_,
        pg.b1s_,
        pg.b2s_,
        pg.b3s_,
        pg.b4s_
    );
}


template<class Specie>
Foam::Tait<Specie> Foam::operator==
(
    const Tait<Specie>& pg1,
    const Tait<Specie>& pg2
)
{
    Specie sp			// unrelevant da keine Mischung
    (
        static_cast<const Specie&>(pg1)
     == static_cast<const Specie&>(pg2)
    );

    const scalar Y1 = pg1.Y()/sp.Y();
    const scalar Y2 = pg2.Y()/sp.Y();


    return Tait<Specie>
    (
        sp
    );
}


// ************************************************************************* //
