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

Application
    rhoSimpleFoam

Description
    Steady-state solver for turbulent flow of compressible fluids.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fluidThermo.H"		//vgl. Erstellung des Felds Thermo in createFields.H	|	Include a thermophysic library
#include "turbulentFluidThermoModel.H"
#include "simpleControl.H"		//Solution control using SIMPLE class
#include "pressureControl.H"		//Provides controls for the pressure reference is closed-volume simulations and a general method for limiting the pressure during the startup of steady-state simulations.
#include "fvOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"		//Create time (object runtime)
    #include "createMesh.H"		//Create time (object mesh)
    #include "createControl.H"
    #include "createFields.H"		//Initialize fields
    #include "createFieldRefs.H"
    #include "CourantNo.H"		//Calculates and outputs Courant Number
    #include "initContinuityErrs.H"	//Declare and initialize the cumulative continuity error

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (simple.loop(runTime))
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;	//Ausgabe welcher Iterationsschritt

        #include "CourantNo.H"		//Calculates and outputs Courant Number

        // Pressure-velocity SIMPLE corrector
        #include "UEqn.H"
        #include "EEqn.H"		//Include an energy conservation equation

        if (simple.consistent())	//vgl. die Abfrage im SIMPLE subdict in fvSolution
        {
            #include "pcEqn.H"		//Verwendung wenn consistent = yes --> d.h. SIMPLEC-Algorithmus
        }
        else
        {
            #include "pEqn.H"		//Verwendung wenn consistent = no --> d.h. SIMPLE-Algorithmus
        }

        turbulence->correct();		// Zugriff auf die virtuelle "correct"-Funktion, Definition in generalizedNewtonian.C	

        runTime.write();
	if(runTime.outputTime())
	{
	//MeinMu.write();	//Tim
	Lambda.write();	//Tim
	alphaAusgabe.write();	//Tim
	}
	//Info<<"_________________________________________________________________________________Aktuelles Kappa = " << thermo.kappa() << nl << endl;
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"		//Ausgabe im Terminal
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"		//Ausgabe im Terminal
            << nl << endl;							//Ausgabe im Terminal

    }

    Info<< "End\n" << endl;		//Ausgabe im Terminal

    return 0;
}


// ************************************************************************* //
