volScalarField muEff = rho*(turbulence->nu() + turbulence->nut());
muEff *= voidfraction;

fvVectorMatrix UEqn
(
    fvm::ddt(voidfractionRho, U) 
  - fvm::Sp(fvc::ddt(voidfractionRho) + fvc::div(rhoPhi), U)
  + fvm::div(rhoPhi, U)
  - fvm::laplacian(muEff, U) - fvc::div(muEff*dev2(fvc::grad(U)().T()))
  ==
    fvOptions(rho,U)
  + rho*g
);

UEqn.relax();
fvOptions.constrain(UEqn);

if (piso.momentumPredictor())
{
	solve
	(
		UEqn
	 ==
		fvc::reconstruct
		(
			(
                 mixture.surfaceTensionForce()
                - voidfractionf*fvc::snGrad(p)
			) * mesh.magSf()
		)
	);
	fvOptions.correct(U);
}