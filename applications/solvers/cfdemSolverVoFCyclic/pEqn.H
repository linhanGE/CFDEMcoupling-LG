{
	volScalarField rAUc(1.0/UEqn.A());
	surfaceScalarField rAUcf(fvc::interpolate(rAUc));

	surfaceScalarField phif
	(
		fvc::flux(rAUc*f)
	);
	
	volVectorField HbyA(constrainHbyA(rAUc*UEqn.H(), U, p));

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::flux(HbyA)
      + voidfractionf*fvc::interpolate(rho*rAUc)*fvc::ddtCorr(U, phi)
    );
	
	adjustPhi(phiHbyA, U, p);
	
    surfaceScalarField phiSt
    (
       mixture.surfaceTensionForce()*rAUcf*mesh.magSf()
    );
	
    phiHbyA += phiSt + phif;
	
	surfaceScalarField rAUcfA = rAUcf*voidfractionf;

	constrainPressure(p, U, phiHbyA, rAUcfA);
	
    while (piso.correctNonOrthogonal())
    {
        surfaceScalarField Dp("Dp", voidfractionf*rAUcfA);

        fvScalarMatrix pEqn
        (
			fvm::laplacian(Dp, p)
         ==
            fvc::ddt(voidfraction) + fvc::div(voidfractionf*phiHbyA)
        );

		pEqn.setReference(pRefCell, getRefCellValue(p, pRefCell));

		pEqn.solve(mesh.solver(p.select(piso.finalInnerIter())));

        if (piso.finalNonOrthogonalIter())
        {
			phi = phiHbyA - pEqn.flux()/voidfractionf;
			
			U = 
			    HbyA 
			  + voidfraction*rAUc*fvc::reconstruct((-pEqn.flux()/voidfractionf)/rAUcfA)
			  + rAUc*fvc::reconstruct(phiSt/rAUcfA)
			  + rAUc*f;

			U.correctBoundaryConditions();
			
            fvOptions.correct(U);
        }
    }
	
	voidfractionPhic = voidfractionf*phi;

    #include "continuityErrs.H"

    if (p.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
    }
}
