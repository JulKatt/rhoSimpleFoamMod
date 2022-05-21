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
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
Foam::Tait<Specie>::Tait
(
    const dictionary& dict
)
//hier u.a. Definition woher welche Parameter eingelesen werden, und auf welche Variable sie intern geschrieben werden
//Aufgrund der Struktur von C wird dies auch in der x.H und xI.H Datei wiederholt, d.h. neue Variablen undbedingt an allen Stellen ergänzen 
//Bsp: b5_(readScalar(dict.subDict("equationOfState").lookup("b5"))),
//	b5_ = interner Variablenname
//	readScalar = Funktion
//	dict.subDict("equationOfState") = diesen Parameter im subDict "equationOfState" suchen (in Datei thermophysicalProperties)
//	lookup("b5") = in dem Subdict "equationOfState" wird nach der Variable b5 gesucht
:
    Specie(dict),
    b5_(readScalar(dict.subDict("equationOfState").lookup("b5"))),
    b6_(readScalar(dict.subDict("equationOfState").lookup("b6"))),
    b1m_(readScalar(dict.subDict("equationOfState").lookup("b1m"))),
    b2m_(readScalar(dict.subDict("equationOfState").lookup("b2m"))),
    b3m_(readScalar(dict.subDict("equationOfState").lookup("b3m"))),
    b4m_(readScalar(dict.subDict("equationOfState").lookup("b4m"))),
    b1s_(readScalar(dict.subDict("equationOfState").lookup("b1s"))),
    b2s_(readScalar(dict.subDict("equationOfState").lookup("b2s"))),
    b3s_(readScalar(dict.subDict("equationOfState").lookup("b3s"))),
    b4s_(readScalar(dict.subDict("equationOfState").lookup("b4s")))
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class Specie>
void Foam::Tait<Specie>::write(Ostream& os) const
{
    Specie::write(os);
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class Specie>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const Tait<Specie>& pg
)
{
    pg.write(os);
    return os;
}


// ************************************************************************* //
