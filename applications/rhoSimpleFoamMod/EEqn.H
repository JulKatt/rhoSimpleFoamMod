{ //vgl. hierzu thermophysicalModels-OF-2.3-or-2.4.x.pdf
    volScalarField& he = thermo.he();		//"he" is either the specific internal energy or the specific enthalpy (choice done when preparing a case, in the dictionary constant/thermophysicalProperties)

    fvScalarMatrix EEqn
    (
        fvm::div(phi, he)				//Konvektionsterm
      + (									//jetzt folgt eine IF-THEN-ELSE Bedingung
            he.name() == "e"							//IF		welches "Modell"?		e = spezifische innere Energie		h = spezifische Enthalpie
          ? fvc::div(phi, volScalarField("Ekp", 0.5*magSqr(U) + p/rho))		//THEN		Für Verwendung innere Energie 		(he = "e")	d.h. Verwendung von sensibleinternalENERGY 
          : fvc::div(phi, volScalarField("K", 0.5*magSqr(U)))			//ELSE		Für Verwendung spezifische Enthalpie 	(he = "h")
        )
      - fvm::laplacian(turbulence->alphaEff(), he)	//alphaEff = alpha laminar (aus src/thermophysicalModels) + alpha turbulent (aus src/turbulenceModels)
     ==
        fvOptions(rho, he)				//Zusätzliche Quell-/Senkterme oder Zwangsbedingungen aus der optionalen Datei fvOptions in case/constant/ 
    );							//Beim starten des Case: "No finite volume options present"  -->  keine Verwendung von fvOptions
//_____________________________________________________________

    EEqn.relax();			//relaxieren 

    fvOptions.constrain(EEqn);		//entfällt da kein fvOptions

    EEqn.solve();			//lösen

    fvOptions.correct(he);		//entfällt da kein fvOptions

    thermo.correct();			// in heRhoThermo, dort als Wrapper --> Aufruf von calculate()
}


//Wenn in obiger Gleichung he = e, wird durch die IF-Bedingung der Term "p/rho" zu der Energiegleichung addiert, was wiederrum h = e + p/rho liefert. Dies wäre bei he = h direkt impliziert.
// --> in jedem Fall also ist der 1. Hauptsatz der Thermodynamik Formelseitig erfüllt.
