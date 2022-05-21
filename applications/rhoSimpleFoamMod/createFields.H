Info<< "Reading thermophysical properties\n" << endl;

autoPtr<fluidThermo> pThermo		//create the object "thermo"
(					//..
    fluidThermo::New(mesh)		//..	autoPtr<fluidThermo> pThermo Map from cell based data back to patch based data
);					//..
fluidThermo& thermo = pThermo();	//create the object "thermo"


thermo.validate(args.executable(), "h", "e");

volScalarField& p = thermo.p();

//_____________________________________________________________
volScalarField rho			//create scalar field rho
(
    IOobject				//create object for input/output operations
    (
        "rho",				//name of the dictionary file to read/write
        runTime.timeName(),		//runtime directory
        mesh,				//object registry
        IOobject::READ_IF_PRESENT,	//Read the dictionary in the runtime directory IF PRESENT
        IOobject::AUTO_WRITE		//Write the value in the runtime directory
    ),
    //mesh				//Link object to mesh
    thermo.rho()			//Quelle
);
//_____________________________________________________________
/*volScalarField mu
(
    IOobject
    (
        "mu",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    //mesh
    mixture.mu()//thermo.mu()
);*/
//_____________________________________________________________
Info<< "Reading field U\n" << endl;
volVectorField U			//create vector field U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh				//Link object to mesh
);

//_____________________________________________________________
#include "compressibleCreatePhi.H"				//Beschreibung für "createPhi.H": creates and initializes the relative face-flux field phi

pressureControl pressureControl(p, rho, simple.dict());

mesh.setFluxRequired(p.name());
//_____________________________________________________________
Info<< "Creating turbulence model\n" << endl;
autoPtr<compressible::turbulenceModel> turbulence		//turbulence model is created  --> Foam::compressible::turbulenceModel  steckt in turbulentFluidThermoModel.H
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);
//_____________________________________________________________

//_____________________________________________________________Tim: Ausgabe des Feldes Kappa (Wärmeleitfähigkeit)
volScalarField Lambda = thermo.kappa();
//Lambda.write();
volScalarField alphaAusgabe = thermo.alpha();
//volScalarField MeinMu = thermo.mu();
//MeinMu.write();

//volScalarField alphah = thermo.alphah();
//alphah.write();

//_____________________________________________________________Ende Tim

dimensionedScalar initialMass = fvc::domainIntegrate(rho);

#include "createMRF.H"
#include "createFvOptions.H"
