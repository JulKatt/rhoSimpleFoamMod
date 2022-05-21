/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::hTaitThermo<EquationOfState>::hTaitThermo
(
    const EquationOfState& st,
    const scalar cp
)
:
    EquationOfState(st),
    Cp_(cp)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::hTaitThermo<EquationOfState>::hTaitThermo
(
    const word& name,
    const hTaitThermo& ct
)
:
    EquationOfState(name, ct),
    Cp_(ct.Cp_)
{}


template<class EquationOfState>
inline Foam::autoPtr<Foam::hTaitThermo<EquationOfState>>
Foam::hTaitThermo<EquationOfState>::clone() const
{
    return autoPtr<hTaitThermo<EquationOfState>>
    (
        new hTaitThermo<EquationOfState>(*this)
    );
}


template<class EquationOfState>
inline Foam::autoPtr<Foam::hTaitThermo<EquationOfState>>
Foam::hTaitThermo<EquationOfState>::New(const dictionary& dict)
{
    return autoPtr<hTaitThermo<EquationOfState>>
    (
        new hTaitThermo<EquationOfState>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::scalar Foam::hTaitThermo<EquationOfState>::limit
(
    const scalar T
) const
{
    return T;
}


template<class EquationOfState>
inline Foam::scalar Foam::hTaitThermo<EquationOfState>::Cp
(
    const scalar p,
    const scalar T
) const
{

    int i_index;
    scalar dT=1;
    scalar T0=293.15;
    scalar Temp_PS;
    scalar Cp_T_PS=0;

    //Überprüfung ob T innerhalb des sinvollen Intervalls (abgesichert durch fvOptions)
    if (T>700 || T<273.15)
        {
            FatalErrorIn
            (
                " Attempt to use Cp-Interpolation,"
                " out of temperature range "
            )   << " T_Cp="<<T<<"; "
                << abort(FatalError);
        }

     // Spezifische Wärmekapazität bei konstantem Druck Cp [J/(kg*K)]
     // Relevant für die Berechnung von T nach dh/dT=Cp 

       scalar Cp_PS[381]=	//Array mit allen 221 Messwerten im Intervall 20°C bis 400°C
     {
	1040,
	1045,
	1050,
	1055,
	1060,
	1065,
	1070,
	1070,
	1070,
	1075,
	1080,
	1085,
	1090,
	1095,
	1100,
	1105,
	1110,
	1110,
	1110,
	1115,
	1120,
	1125,
	1130,
	1135,
	1140,
	1145,
	1150,
	1150,
	1150,
	1155,
	1160,
	1165,
	1170,
	1175,
	1180,
	1185,
	1190,
	1195,
	1200,
	1205,
	1210,
	1210,
	1210,
	1215,
	1220,
	1225,
	1230,
	1235,
	1240,
	1245,
	1250,
	1255,
	1260,
	1260,
	1260,
	1265,
	1270,
	1275,
	1280,
	1280,
	1280,
	1285,
	1290,
	1290,
	1290,
	1295,
	1300,
	1305,
	1310,
	1315,
	1320,
	1320,
	1320,
	1325,
	1330,
	1335,
	1340,
	1340,
	1340,
	1345,
	1350,
	1355,
	1360,
	1365,
	1370,
	1375,
	1380,
	1385,
	1390,
	1405,
	1420,
	1420,
	1420,
	1425,
	1430,
	1435,
	1440,
	1450,
	1460,
	1460,
	1460,
	1465,
	1470,
	1475,
	1480,
	1480,
	1480,
	1485,
	1490,
	1490,
	1490,
	1495,
	1500,
	1500,
	1500,
	1505,
	1510,
	1510,
	1510,
	1515,
	1520,
	1525,
	1530,
	1530,
	1530,
	1530,
	1530,
	1535,
	1540,
	1545,
	1550,
	1550,
	1550,
	1555,
	1560,
	1565,
	1570,
	1570,
	1570,
	1575,
	1580,
	1580,
	1580,
	1585,
	1590,
	1590,
	1590,
	1595,
	1600,
	1600,
	1600,
	1605,
	1610,
	1615,
	1620,
	1615,
	1610,
	1620,
	1630,
	1625,
	1620,
	1630,
	1640,
	1645,
	1650,
	1655,
	1660,
	1660,
	1660,
	1665,
	1670,
	1670,
	1670,
	1675,
	1680,
	1680,
	1680,
	1685,
	1690,
	1695,
	1700,
	1700,
	1700,
	1705,
	1710,
	1710,
	1710,
	1715,
	1720,
	1725,
	1730,
	1730,
	1730,
	1735,
	1740,
	1740,
	1740,
	1745,
	1750,
	1750,
	1750,
	1755,
	1760,
	1760,
	1760,
	1765,
	1770,
	1770,
	1770,
	1775,
	1780,
	1780,
	1780,
	1785,
	1790,
	1790,
	1790,
	1790,
	1790,
	1790,
	1790,
	1790,
	1790,
	1790,
	1790,
	1790,
	1790,
	1790,
	1790,
	1795,
	1800,
	1800,
	1800,
	1800,
	1800,
	1805,
	1810,
	1810,
	1810,
	1810,
	1810,
	1815,
	1820,
	1820,
	1820,
	1820,
	1820,
	1825,
	1830,
	1830,
	1830,
	1830,
	1830,
	1835,
	1840,
	1840,
	1840,
	1845,
	1850,
	1850,
	1850,
	1850,
	1850,
	1855,
	1860,
	1860,
	1860,
	1865,
	1870,
	1870,
	1870,
	1875,
	1880,
	1885,
	1890,
	1890,
	1890,
	1895,
	1900,
	1900,
	1900,
	1905,
	1910,
	1915,
	1920,
	1925,
	1930,
	1930,
	1930,
	1935,
	1940,
	1945,
	1950,
	1955,
	1960,
	1970,
	1980,
	1985,
	1990,
	2005,
	2020,
	2040,
	2060,
	2140,
	2220,
	2525,
	2830,
	2970,
	3110,
	3060,
	3010,
	2810,
	2610,
	2380,
	2150,
	2075,
	2000,
	1995,
	1990,
	1985,
	1980,
	1985,
	1990,
	1990,
	1990,
	1990,
	1990,
	1990,
	1990,
	1995,
	2000,
	2000,
	2000,
	2005,
	2010,
	2010,
	2010,
	2010,
	2010,
	2015,
	2020,
	2025,
	2030,
	2030,
	2030,
	2030,
	2030,
	2035,
	2040,
	2045,
	2050,
	2055,
	2060,
	2060,
	2060,
	2065,
	2070,
	2075,
	2080,
	2085,
	2090,
	2090,
	2090,
	2095,
	2100,
	2105,
	2110,
	2115,
	2120,
	2120,
	2120,
	2125,
	2130,
	2135,
	2140,
	2140,
	2140,
	2140,
	2140,
	2145,
	2150
	};
         
     //_______________________________lineare interpolation für Cp(T)

    if (T<293.15)	//konstanter Wert für alles kleiner als dem Messwert bei der kleinsten Temperatur
        {
   		Cp_T_PS = 1040;
        }

    if (T>=293.15&& T<673.15)	//lineare Interpolation im Messbereich 20°C bis 240°C
        {
     		i_index = int(floor(fabs((T-T0)/dT)));
     		Temp_PS = T0+i_index*dT;
   		Cp_T_PS = Cp_PS[i_index] + (Cp_PS[i_index+1]-Cp_PS[i_index])*(T-Temp_PS)/dT;
        }

    if (T>=673.15)	//konstanter Wert für alles größer als dem Messwert bei der höchsten Temperatur
        {
   		Cp_T_PS = 2150;
        }

    if (Cp_T_PS>3110|| Cp_T_PS<1040)	//Überprüfung ob Werte außerhalb des definierten Intervalls
        {
		Info<< "Fehler bei Cp" << Cp_T_PS << endl;
        }

     return Cp_T_PS;	//Cp(T) könnte auch in der EOS TaitI.H implementiert werden, dafür hier "return EquationOfState::Cp(p, T);" um darauf zuzugreifen
}


template<class EquationOfState>
inline Foam::scalar Foam::hTaitThermo<EquationOfState>::Hs
(
    const scalar p, const scalar T
) const
{

    int i_index;
    scalar dT=1;
    scalar T0=293.15;
    scalar Temp_PS;
    scalar Hs_T_PS=0;

    //Überprüfung ob T innerhalb des sinvollen Intervalls (abgesichert durch fvOptions)
    if (T>700 || T<273.15)
        {
            FatalErrorIn
            (
                " Attempt to use Hs-Interpolation,"
                " out of temperature range "
            )   << " T_Hs="<<T<<"; "
                << abort(FatalError);
        }

      scalar Hs_PS[381]=	//Array mit allen 221 Messwerten im Intervall 20°C bis 240°C
     {
	0,
	1045,
	2090,
	3145,
	4200,
	5265,
	6330,
	7400,
	8470,
	9545,
	10620,
	11705,
	12790,
	13885,
	14980,
	16085,
	17190,
	18300,
	19410,
	20525,
	21640,
	22765,
	23890,
	25025,
	26160,
	27305,
	28450,
	29600,
	30750,
	31905,
	33060,
	34225,
	35390,
	36565,
	37740,
	38925,
	40110,
	41305,
	42500,
	43705,
	44910,
	46120,
	47330,
	48545,
	49760,
	50985,
	52210,
	53445,
	54680,
	55925,
	57170,
	58425,
	59680,
	60940,
	62200,
	63465,
	64730,
	66005,
	67280,
	68560,
	69840,
	71125,
	72410,
	73700,
	74990,
	76285,
	77580,
	78885,
	80190,
	81505,
	82820,
	84140,
	85460,
	86785,
	88110,
	89445,
	90780,
	92120,
	93460,
	94805,
	96150,
	97505,
	98860,
	100225,
	101590,
	102965,
	104340,
	105725,
	107110,
	108515,
	109920,
	111340,
	112760,
	114185,
	115610,
	117045,
	118480,
	119930,
	121380,
	122840,
	124300,
	125765,
	127230,
	128705,
	130180,
	131660,
	133140,
	134625,
	136110,
	137600,
	139090,
	140585,
	142080,
	143580,
	145080,
	146585,
	148090,
	149600,
	151110,
	152625,
	154140,
	155665,
	157190,
	158720,
	160250,
	161780,
	163310,
	164845,
	166380,
	167925,
	169470,
	171020,
	172570,
	174125,
	175680,
	177245,
	178810,
	180380,
	181950,
	183525,
	185100,
	186680,
	188260,
	189845,
	191430,
	193020,
	194610,
	196205,
	197800,
	199400,
	201000,
	202605,
	204210,
	205825,
	207440,
	209055,
	210670,
	212290,
	213910,
	215535,
	217160,
	218790,
	220420,
	222065,
	223710,
	225365,
	227020,
	228680,
	230340,
	232005,
	233670,
	235340,
	237010,
	238685,
	240360,
	242040,
	243720,
	245405,
	247090,
	248785,
	250480,
	252180,
	253880,
	255585,
	257290,
	259000,
	260710,
	262425,
	264140,
	265865,
	267590,
	269320,
	271050,
	272785,
	274520,
	276260,
	278000,
	279745,
	281490,
	283240,
	284990,
	286745,
	288500,
	290260,
	292020,
	293785,
	295550,
	297320,
	299090,
	300865,
	302640,
	304420,
	306200,
	307985,
	309770,
	311560,
	313350,
	315140,
	316930,
	318720,
	320510,
	322300,
	324090,
	325880,
	327670,
	329460,
	331250,
	333040,
	334830,
	336625,
	338420,
	340220,
	342020,
	343820,
	345620,
	347425,
	349230,
	351040,
	352850,
	354660,
	356470,
	358285,
	360100,
	361920,
	363740,
	365560,
	367380,
	369205,
	371030,
	372860,
	374690,
	376520,
	378350,
	380185,
	382020,
	383860,
	385700,
	387545,
	389390,
	391240,
	393090,
	394940,
	396790,
	398645,
	400500,
	402360,
	404220,
	406085,
	407950,
	409820,
	411690,
	413565,
	415440,
	417325,
	419210,
	421100,
	422990,
	424885,
	426780,
	428680,
	430580,
	432485,
	434390,
	436305,
	438220,
	440145,
	442070,
	444000,
	445930,
	447865,
	449800,
	451745,
	453690,
	455645,
	457600,
	459570,
	461540,
	463525,
	465510,
	467515,
	469520,
	471560,
	473600,
	475740,
	477880,
	480405,
	482930,
	485900,
	488870,
	491930,
	494990,
	497800,
	500610,
	502990,
	505370,
	507445,
	509520,
	511515,
	513510,
	515495,
	517480,
	519465,
	521450,
	523440,
	525430,
	527420,
	529410,
	531400,
	533390,
	535385,
	537380,
	539380,
	541380,
	543385,
	545390,
	547400,
	549410,
	551420,
	553430,
	555445,
	557460,
	559485,
	561510,
	563540,
	565570,
	567600,
	569630,
	571665,
	573700,
	575745,
	577790,
	579845,
	581900,
	583960,
	586020,
	588085,
	590150,
	592225,
	594300,
	596385,
	598470,
	600560,
	602650,
	604745,
	606840,
	608945,
	611050,
	613165,
	615280,
	617400,
	619520,
	621645,
	623770,
	625905,
	628040,
	630180,
	632320,
	634460,
	636600,
	638745,
	640890
    	};
         
     //_______________________________lineare interpolation für Hs(T)

    if (T<293.15)	//konstanter Wert für alles kleiner als dem Messwert bei der kleinsten Temperatur
        {
   		Hs_T_PS = 0;
        }

    if (T>=293.15&& T<673.15)	//lineare Interpolation im Messbereich 20°C bis 240°C
        {
     		i_index = int(floor(fabs((T-T0)/dT)));
     		Temp_PS = T0+i_index*dT;
   		Hs_T_PS = Hs_PS[i_index] + (Hs_PS[i_index+1]-Hs_PS[i_index])*(T-Temp_PS)/dT;
        }

    if (T>=673.15)	//konstanter Wert für alles größer als dem Messwert bei der höchsten Temperatur
        {
   		Hs_T_PS = 640890;
        }

    if (Hs_T_PS>640890|| Hs_T_PS<0)	//Überprüfug ob Wert außerhalb des definierten Intervalls
        {
		Info<< "Fehler bei Hs" << Hs_T_PS << endl;
        }

     return Hs_T_PS;
}



template<class EquationOfState>
inline Foam::scalar Foam::hTaitThermo<EquationOfState>::Hc() const
{
    return 0;		//Heat of Fusion, d.h. Bildungsenthalpie --> nur für Mischungen relevant --> 0
}


template<class EquationOfState>
inline Foam::scalar Foam::hTaitThermo<EquationOfState>::Ha
(
    const scalar p, const scalar T
) const
{
    return Hs(p, T) + Hc();	//Hc()=0 siehe oben
}


template<class EquationOfState>
inline Foam::scalar Foam::hTaitThermo<EquationOfState>::S
(
    const scalar p, const scalar T
) const
{
	//Info<< "_____________________________________________________________________Aufruf" << endl;
    return 0;
}


template<class EquationOfState>
inline Foam::scalar Foam::hTaitThermo<EquationOfState>::dGdT
(
    const scalar p, const scalar T
) const
{
    return 0;
}


template<class EquationOfState>
inline Foam::scalar Foam::hTaitThermo<EquationOfState>::dCpdT
(
    const scalar p, const scalar T
) const
{
    return 0;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class EquationOfState>
inline void Foam::hTaitThermo<EquationOfState>::operator+=
(
    const hTaitThermo<EquationOfState>& ct
)
{
    scalar Y1 = this->Y();

    EquationOfState::operator+=(ct);

    if (mag(this->Y()) > small)
    {
        Y1 /= this->Y();
        scalar Y2 = ct.Y()/this->Y();

        Cp_ = Y1*Cp_ + Y2*ct.Cp_;
    }
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::hTaitThermo<EquationOfState> Foam::operator+
(
    const hTaitThermo<EquationOfState>& ct1,
    const hTaitThermo<EquationOfState>& ct2
)
{
    EquationOfState eofs
    (
        static_cast<const EquationOfState&>(ct1)
      + static_cast<const EquationOfState&>(ct2)
    );

    if (mag(eofs.Y()) < small)
    {
        return hTaitThermo<EquationOfState>
        (
            eofs,
            ct1.Cp_
        );
    }
    else
    {
        return hTaitThermo<EquationOfState>
        (
            eofs,
            ct1.Y()/eofs.Y()*ct1.Cp_
          + ct2.Y()/eofs.Y()*ct2.Cp_
        );
    }
}


template<class EquationOfState>
inline Foam::hTaitThermo<EquationOfState> Foam::operator*
(
    const scalar s,
    const hTaitThermo<EquationOfState>& ct
)
{
    return hTaitThermo<EquationOfState>
    (
        s*static_cast<const EquationOfState&>(ct),
        ct.Cp_
    );
}


template<class EquationOfState>
inline Foam::hTaitThermo<EquationOfState> Foam::operator==
(
    const hTaitThermo<EquationOfState>& ct1,
    const hTaitThermo<EquationOfState>& ct2
)
{
    EquationOfState eofs
    (
        static_cast<const EquationOfState&>(ct1)
     == static_cast<const EquationOfState&>(ct2)
    );

    return hTaitThermo<EquationOfState>
    (
        eofs,
        ct2.Y()/eofs.Y()*ct2.Cp_
      - ct1.Y()/eofs.Y()*ct1.Cp_
    );
}


// ************************************************************************* //
