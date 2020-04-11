#include "createRDeltaT.H"

Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volVectorField Usmooth
(
    IOobject
    (
        "Usmooth",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh.lookupObject<volVectorField> ("U")
);

#include "createPhi.H"

//===============================
// particle interaction modelling
//===============================

Info<< "\nReading momentum exchange field f\n" << endl;
volVectorField f
    (
        IOobject
        (
            "f",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
        //dimensionedScalar("0", dimensionSet(1, -3, -1, 0, 0), 1.0)
    );

volScalarField particlefraction
(
    IOobject
    (
        "particlefraction",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// voidfraction must be constructed before the cloud
// so that the drag-models can find it
volScalarField voidfraction
(
    IOobject
    (
        "voidfraction",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh
);
voidfraction.oldTime();

Info<< "Reading transportProperties\n" << endl;
immiscibleIncompressibleTwoPhaseMixture mixture(U, phi);

volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());

const dimensionedScalar& rho1 = mixture.rho1();
const dimensionedScalar& rho2 = mixture.rho2();

// Need to store rho for ddt(rho, U)
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    alpha1*rho1 + alpha2*rho2
);
rho.oldTime();            // automatically recover the old time rho

Info<< "Reading particle velocity field Us\n" << endl;
volVectorField Us
(
	IOobject
	(
		"Us",
		runTime.timeName(),
		mesh,
		IOobject::MUST_READ,
		IOobject::AUTO_WRITE
	),
	mesh
);

// Need to store mu as incompressibleTwoPhaseMixture does not store it
volScalarField mu
(
    IOobject
    (
        "mu",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT
    ),
    mixture.mu(),
    calculatedFvPatchScalarField::typeName
);

// Mass flux
 surfaceScalarField rhoPhi
 (
     IOobject
     (
         "rhoPhi",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::NO_WRITE
     ),
     fvc::flux(rho*U*voidfraction)
 );

#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"

volScalarField sigmaK
(
    IOobject
    (
        "K",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mixture.sigmaK()
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    mesh.solutionDict().subDict("PISO"),
    pRefCell,
    pRefValue
);

if (p.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
}

mesh.setFluxRequired(p.name());
mesh.setFluxRequired(alpha1.name());

volScalarField voidfractionRho(voidfraction*rho);
voidfractionRho.oldTime();

surfaceScalarField voidfractionf("voidfractionf", fvc::interpolate(voidfraction));

// Volumetric phase flux
surfaceScalarField voidfractionPhic("voidfractionPhic", voidfractionf*phi);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, mixture)
);