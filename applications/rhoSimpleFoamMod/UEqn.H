    // Solve the Momentum equation

    MRF.correctBoundaryVelocity(U);		//keine MRF

    tmp<fvVectorMatrix> tUEqn
    (
        fvm::div(phi, U)			//Konvektionsterm
      + MRF.DDt(rho, U)				//entfällt da keine Verwendung der Multi-Reference-Frame (MRF) Annahme
      + turbulence->divDevRhoReff(U)		//Objekt "turbulence" ruft Funktion "divDevRhoReff" auf mit dem Argument "U"
     ==
        fvOptions(rho, U)			//entfällt da keine zusätzlichen Quell-/Senkterme in optionaler Datei fvOptions
    );
    fvVectorMatrix& UEqn = tUEqn.ref();

    UEqn.relax();				//relaxieren

    fvOptions.constrain(UEqn);			//entfällt

    solve(UEqn == -fvc::grad(p));		//lösen gegen den Druck

    fvOptions.correct(U);			//entfällt
