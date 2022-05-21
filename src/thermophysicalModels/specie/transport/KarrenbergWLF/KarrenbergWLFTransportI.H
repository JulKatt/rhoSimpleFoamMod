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

\*---------------------------------------------------------------------------*/

#include "specie.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::KarrenbergWLFTransport<Thermo>::KarrenbergWLFTransport
(
    const word& name,
    const KarrenbergWLFTransport& Kwlft
)
:
    Thermo(name, Kwlft),
    mu0_(Kwlft.mu0_),
    Tr_(Kwlft.Tr_),
    Tm_(Kwlft.Tm_),
    TS_(Kwlft.TS_),
    MS_(Kwlft.MS_),
    C1_(Kwlft.C1_),
    C2_(Kwlft.C2_)
{}


template<class Thermo>
inline Foam::autoPtr<Foam::KarrenbergWLFTransport<Thermo>>
Foam::KarrenbergWLFTransport<Thermo>::clone() const
{
    return autoPtr<KarrenbergWLFTransport<Thermo>>
    (
        new KarrenbergWLFTransport<Thermo>(*this)
    );
}


template<class Thermo>
inline Foam::autoPtr<Foam::KarrenbergWLFTransport<Thermo>>
Foam::KarrenbergWLFTransport<Thermo>::New
(
    const dictionary& dict
)
{
    return autoPtr<KarrenbergWLFTransport<Thermo>>
    (
        new KarrenbergWLFTransport<Thermo>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::scalar Foam::KarrenbergWLFTransport<Thermo>::mu
(
    const scalar p,
    const scalar T
) const
{
	const scalar Tg = 373.15;	//Tg ist hier zusätzlicher Parameter für die Modellierung, legt fest ab wo mit konstanter Nullviskosität gerechnet wird

	if (((-C1_*(Tm_ - Tr_)/(C2_ + Tm_ - Tr_))) < (0))
	{
		Info << "Warnung!! Anomalie im WLF-Modell" << endl;
	}


	if ((T-Tg) < (0))
	{
		return mu0_*exp((1-this->MF(p, T))*(((-C1_*(Tm_ - Tr_)/(C2_ + Tm_ - Tr_))))+this->MF(p, T)*((-C1_*(Tg - Tr_)/(C2_ + Tg - Tr_))))+10e8*(-this->MF(p, T)+1);	//Übergang von Tm zu Tg
	}
	else
	{
		return mu0_*exp((1-this->MF(p, T))*(((-C1_*(Tm_ - Tr_)/(C2_ + Tm_ - Tr_))))+this->MF(p, T)*((-C1_*(T - Tr_)/(C2_ + T - Tr_))))+10e8*(-this->MF(p, T)+1);		//Übergang von Tg zu T (basierend auf Tm)
	}

}

template<class Thermo>
inline Foam::scalar Foam::KarrenbergWLFTransport<Thermo>::MF	////_______Berechnung der Melt-Fraction nach dem Karrenberg-Ansatz
(
    const scalar p, const scalar T
) const
{
	return Foam::pow(((Foam::tanh((T-Tm_)*TS_))+1.0)/2.0,MS_);
}


template<class Thermo>
inline Foam::scalar Foam::KarrenbergWLFTransport<Thermo>::kappa
(
    const scalar p, const scalar T
) const
{
    //Info << "Hier Kappa" << endl;			//	--> wird nicht aufgerufen
    return this->Cp(p, T)*mu(p, T);
}


template<class Thermo>
inline Foam::scalar Foam::KarrenbergWLFTransport<Thermo>::alphah	//alphah = thermal diffusivity of enthalpy  --> alphah = (kappa/Cp)	https://openfoam.top/en/alpha/
(
    const scalar p,
    const scalar T
) const
{		//alphah = temperaturleitfähigkeit*Dichte

    int i_index = 0;
    scalar dT=1;
    scalar T0=296.15;		//Messwerte hier erst ab 23°C
    scalar Temp_PS = 0;
    scalar alpha_T_PS=0;

    //Überprüfung ob T innerhalb des sinvollen Intervalls (abgesichert durch fvOptions)
    if (T>700 || T<273.15)
        {
            FatalErrorIn
            (
                " Attempt to use alpha-Interpolation,"
                " out of temperature range "
            )   << " T_alpha="<<T<<"; "
                << abort(FatalError);
        }

      scalar alpha_PS[218]=	//Array mit allen 221 Messwerten im Intervall 20°C bis 240°C
     {
	0.12367,	//Index:0	Temp:23°C
	0.12352,
	0.12338,
	0.12323,
	0.1231,
	0.12296,
	0.12283,
	0.1227,
	0.12257,
	0.12244,
	0.12232,
	0.1222,
	0.12208,
	0.12197,
	0.12185,
	0.12174,
	0.12164,
	0.12153,
	0.12143,
	0.12133,
	0.12124,
	0.12114,
	0.12105,
	0.12096,
	0.12088,
	0.12079,
	0.12071,
	0.12063,
	0.12056,
	0.12048,
	0.12041,
	0.12035,
	0.12028,
	0.12022,
	0.12016,
	0.1201,
	0.12005,
	0.12,
	0.11995,
	0.1199,
	0.11986,
	0.11981,
	0.11977,
	0.11973,
	0.1197,
	0.11967,
	0.11963,
	0.11961,
	0.11958,
	0.11956,
	0.11954,
	0.11952,
	0.11951,
	0.1195,
	0.11949,
	0.11948,
	0.11948,
	0.11948,
	0.11948,
	0.11949,
	0.1195,
	0.11951,
	0.11953,
	0.11955,
	0.11958,
	0.1196,
	0.11963,
	0.11967,
	0.12003,
	0.12092,
	0.12219,
	0.12373,
	0.12538,
	0.12702,
	0.1285,
	0.12969,
	0.13046,
	0.13066,
	0.13034,
	0.12966,
	0.12871,
	0.12756,
	0.12632,
	0.12506,
	0.12388,
	0.12286,
	0.12209,
	0.12166,
	0.12141,
	0.12116,
	0.12092,
	0.12069,
	0.12049,
	0.12034,
	0.12023,
	0.12019,
	0.12022,
	0.12034,
	0.1205,
	0.12067,
	0.12084,
	0.12103,
	0.12122,
	0.12142,
	0.12163,
	0.12185,
	0.12209,
	0.12234,
	0.12261,
	0.12291,
	0.12324,
	0.12357,
	0.12391,
	0.12424,
	0.12455,
	0.12485,
	0.12511,
	0.12534,
	0.12553,
	0.12572,
	0.12589,
	0.12605,
	0.1262,
	0.12634,
	0.12647,
	0.12659,
	0.12671,
	0.12681,
	0.1269,
	0.12699,
	0.12706,
	0.12713,
	0.12719,
	0.12725,
	0.1273,
	0.12734,
	0.12738,
	0.12741,
	0.12743,
	0.12745,
	0.12747,
	0.12748,
	0.12748,
	0.12749,
	0.12749,
	0.12748,
	0.12748,
	0.12747,
	0.12746,
	0.12745,
	0.12743,
	0.12742,
	0.1274,
	0.12739,
	0.12737,
	0.12736,
	0.12735,
	0.12733,
	0.12732,
	0.1273,
	0.12728,
	0.12726,
	0.12723,
	0.1272,
	0.12717,
	0.12714,
	0.1271,
	0.12706,
	0.12701,
	0.12696,
	0.12691,
	0.12686,
	0.1268,
	0.12674,
	0.12668,
	0.12661,
	0.12654,
	0.12647,
	0.12639,
	0.12631,
	0.12623,
	0.12614,
	0.12606,
	0.12596,
	0.12587,
	0.12577,
	0.12567,
	0.12557,
	0.12546,
	0.12535,
	0.12524,
	0.12512,
	0.125,
	0.12488,
	0.12475,
	0.12462,
	0.12449,
	0.12436,
	0.12422,
	0.12408,
	0.12393,
	0.12378,
	0.12363,
	0.12348,
	0.12332,
	0.12316,
	0.123,
	0.12283,
	0.12266,
	0.12249,
	0.12232,
	0.12214,
	0.12196,
	0.12177,
	0.12158,
	0.12139,
	0.1212,
	0.121
     };
         
     //_______________________________lineare interpolation für alpha(T)

    if (T<296.15)	//konstanter Wert für alles kleiner als dem Messwert bei der kleinsten Temperatur
        {
   		alpha_T_PS = 0.12367;
        }

    if (T>=296.15&& T<513.15)	//lineare Interpolation im Messbereich 20°C bis 240°C
        {
     		//i_index = int(floor(fabs((T-T0)/dT)));
     		i_index = int(floor((T-T0)/dT));
     		Temp_PS = T0+i_index*dT;
   		alpha_T_PS = alpha_PS[i_index] + (alpha_PS[i_index+1]-alpha_PS[i_index])*(T-Temp_PS)/dT;
        }

    if (T>=513.15)	//konstanter Wert für alles größer als dem Messwert bei der höchsten Temperatur
        {
   		alpha_T_PS = 0.2;	//0.121;
        }

    if (alpha_T_PS<0.11948 || alpha_T_PS>0.3)	//Überprüfug ob Wert außerhalb des definierten Intervalls, mit Debugging Ausgaben
        {
	Info<< "Fehler bei alpha " << " alpha_PS[i] " << alpha_PS[i_index] << " alpha_PS[i+1] " << alpha_PS[i_index+1] << " T " << T << " Temp_PS " << Temp_PS << " dT " << dT << alpha_T_PS << endl;
        }
	
     return (alpha_T_PS*1e-6)*(this->rho(p, T));	//Umrechung von Temperaturleitfähigkeit in die richtige Einheit,, anschließend Multiplikation mit rho um alphah zurückzugeben
}							//alphah = thermal diffusivity of enthalphy


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline void Foam::KarrenbergWLFTransport<Thermo>::operator+=
(
    const KarrenbergWLFTransport<Thermo>& Kwlft
)
{
    scalar Y1 = this->Y();

    Thermo::operator+=(Kwlft);

    if (mag(this->Y()) > small)
    {
        Y1 /= this->Y();
        scalar Y2 = Kwlft.Y()/this->Y();

        mu0_ = Y1*mu0_ + Y2*Kwlft.mu0_;
        Tr_ = Y1*Tr_ + Y2*Kwlft.Tr_;
        Tm_ = Y1*Tm_ + Y2*Kwlft.Tm_;
        TS_ = Y1*TS_ + Y2*Kwlft.TS_;
        MS_ = Y1*MS_ + Y2*Kwlft.MS_;
        C1_ = Y1*C1_ + Y2*Kwlft.C1_;
        C2_ = Y1*C2_ + Y2*Kwlft.C2_;
    }
}


template<class Thermo>
inline void Foam::KarrenbergWLFTransport<Thermo>::operator*=
(
    const scalar s
)
{
    Thermo::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::KarrenbergWLFTransport<Thermo> Foam::operator+
(
    const KarrenbergWLFTransport<Thermo>& Kwlft1,
    const KarrenbergWLFTransport<Thermo>& Kwlft2
)
{
    Thermo t
    (
        static_cast<const Thermo&>(Kwlft1) + static_cast<const Thermo&>(Kwlft2)
    );

    if (mag(t.Y()) < small)
    {
        return KarrenbergWLFTransport<Thermo>
        (
            t,
            0,
            Kwlft1.mu0_,
            Kwlft1.Tr_,
            Kwlft1.C1_,
            Kwlft1.C2_
        );
    }
    else
    {
        scalar Y1 = Kwlft1.Y()/t.Y();
        scalar Y2 = Kwlft2.Y()/t.Y();
        return KarrenbergWLFTransport<Thermo>
        (
            t,
            Y1*Kwlft1.mu0_ + Y2*Kwlft2.mu0_,
            Y1*Kwlft1.Tr_ + Y2*Kwlft2.Tr_,
            Y1*Kwlft1.C1_ + Y2*Kwlft2.C1_,
            Y1*Kwlft1.C2_ + Y2*Kwlft2.C2_
        );
    }
}


template<class Thermo>
inline Foam::KarrenbergWLFTransport<Thermo> Foam::operator*
(
    const scalar s,
    const KarrenbergWLFTransport<Thermo>& Kwlft
)
{
    return KarrenbergWLFTransport<Thermo>
    (
        s*static_cast<const Thermo&>(Kwlft),
        Kwlft.mu0_,
        Kwlft.Tr_,
        Kwlft.Tm_,
        Kwlft.TS_,
        Kwlft.MS_,
        Kwlft.C1_,
        Kwlft.C2_
    );
}


// ************************************************************************* //
