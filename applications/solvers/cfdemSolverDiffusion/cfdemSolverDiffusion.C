/*---------------------------------------------------------------------------*\
    CFDEMcoupling - Open Source CFD-DEM coupling

    CFDEMcoupling is part of the CFDEMproject
    www.cfdem.com
                                Christoph Goniva, christoph.goniva@cfdem.com
                                Copyright (C) 1991-2009 OpenCFD Ltd.
                                Copyright (C) 2009-2012 JKU, Linz
                                Copyright (C) 2012-     DCS Computing GmbH,Linz
-------------------------------------------------------------------------------
License
    This file is part of CFDEMcoupling.

    CFDEMcoupling is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    CFDEMcoupling is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with CFDEMcoupling.  If not, see <http://www.gnu.org/licenses/>.

Application
    cfdemSolverDiffusion

Description
    Transient solver for incompressible flow.
    Turbulence modelling is generic, i.e. laminar, RAS or LES may be selected.
    The code is an evolution of the solver pisoFoam in OpenFOAM(R) 1.6,
    where additional functionality for CFD-DEM coupling is added.
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"

#include "OFversion.H"
#if defined(version30)
    #include "turbulentTransportModel.H"
    #include "pisoControl.H"
#else
    #include "turbulenceModel.H"
#endif
#if defined(versionv1606plus) || defined(version40)
    #include "fvOptions.H"
#else
    #include "fvIOoptionList.H"
#endif
#include "fixedFluxPressureFvPatchScalarField.H"
#include "cfdemCloud.H"

#if defined(anisotropicRotation)
    #include "cfdemCloudRotation.H"
#endif
#if defined(superquadrics_flag)
    #include "cfdemCloudRotationSuperquadric.H"
#endif
#include "explicitCouple.H"
#include "clockModel.H"
#include "smoothingModel.H"
#include "forceModel.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #if defined(version30)
        pisoControl piso(mesh);
        #include "createTimeControls.H"
    #endif
    #include "createFields.H"
    #include "createFvOptions.H"
    #include "initContinuityErrs.H"

    // create cfdemCloud
    #include "readGravitationalAcceleration.H"
    // #include "checkImCoupleM.H"
    #if defined(anisotropicRotation)
        cfdemCloudRotation particleCloud(mesh);
    #elif defined(superquadrics_flag)
        cfdemCloudRotationSuperquadric particleCloud(mesh);
    #else
        cfdemCloud particleCloud(mesh);
    #endif
    #include "checkModelType.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info<< "\nStarting time loop\n" << endl;
    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

		#include "readTimeControls.H"
		#include "CourantNo.H"
		#include "setDeltaT.H"

        // do particle stuff
        particleCloud.clockM().start(1,"Global");
        particleCloud.clockM().start(2,"Coupling");
        bool hasEvolved = particleCloud.diffusionEvolve(voidfraction,Us,U,Usmooth);
		
		if(hasEvolved)
		{
            particleCloud.smoothingM().smoothenAbsolutField(particleCloud.forceM(0).expParticleForces());
		}
		
		// be careful with the negative sign, becasue we need the force on fluid
        f = -1 * particleCloud.momCoupleM(particleCloud.registryM().getProperty("explicitCouple_index")).expMomSource();  
        f.correctBoundaryConditions();
		        	
        surfaceScalarField voidfractionf = fvc::interpolate(voidfraction);
        phi = voidfractionf*phiByVoidfraction;

        //Force Checks
        #include "forceCheckEx.H"           // note!

        // #include "solverDebugInfo.H"
        particleCloud.clockM().stop("Coupling");

        particleCloud.clockM().start(26,"Flow");

        if(particleCloud.solveFlow())
        {
            // Pressure-velocity PISO corrector
            {
                // Momentum predictor
                fvVectorMatrix UEqn
                (
                    fvm::ddt(voidfraction,U) - fvm::Sp(fvc::ddt(voidfraction),U)
                  + fvm::div(phi,U) - fvm::Sp(fvc::div(phi),U)
                  + particleCloud.divVoidfractionTau(U, voidfraction)
                  ==
                  fvOptions(U)
                );

                UEqn.relax();
                fvOptions.constrain(UEqn);
               
                // the momentum prediction is not neccessary
                if (piso.momentumPredictor())
                {
                    if (modelType=="B" || modelType=="Bfull")
                        solve(UEqn == - fvc::grad(p));
                    else
                        solve(UEqn == - voidfraction*fvc::grad(p));

                    fvOptions.correct(U);
                }
      
                // --- PISO loop
                while (piso.correct())
                {
                    volScalarField rUA = 1.0/UEqn.A();   

                    surfaceScalarField rUAf("(1|A(U))", fvc::interpolate(rUA));
                    surfaceScalarField rUAfvoidfraction("(voidfraction2|A(U)F)", rUAf*voidfractionf);

                    U = rUA*UEqn.H();    // U = HbyA

                    surfaceScalarField phiForces
					(
						fvc::flux(rUA*f/rho)
					);

                    /*forAll(p.boundaryFieldRef(), patchi)     // why we need this
					{
						if (isA<zeroGradientFvPatchScalarField>(p.boundaryFieldRef()[patchi]))
						{
							phiForces.boundaryFieldRef()[patchi] = 0.0;
						}
					}
					Us.correctBoundaryConditions();*/

                    phi = fvc::flux(U)   // use fvc::flux can reduce the storage use of interpolated field
                          + rUAfvoidfraction*fvc::ddtCorr(U, phiByVoidfraction);
						  
                    phi += phiForces;

                    if (modelType=="A")
                        rUAfvoidfraction = surfaceScalarField("(voidfraction2|A(U))",rUAf*voidfractionf*voidfractionf);

                    // Update the fixedFluxPressure BCs to ensure flux consistency
                    #include "fixedFluxPressureHandling.H"    
					
					// surfaceScalarField phiSF = alphasf*fvc::flux(Us) + voidfractionf*phi;  // total flux of fluid-particle mixture

                    // Non-orthogonal pressure corrector loop
                    while (piso.correctNonOrthogonal())
                    {
                        // Pressure corrector
                        fvScalarMatrix pEqn
                        (
                            fvm::laplacian(rUAfvoidfraction, p) == fvc::div(voidfractionf*phi) + particleCloud.ddtVoidfraction()
                        );

                        pEqn.setReference(pRefCell, pRefValue);

                        pEqn.solve(mesh.solver(p.select(piso.finalInnerIter())));

                        if (piso.finalNonOrthogonalIter())
                        {
                            phiByVoidfraction = phi - pEqn.flux()/voidfractionf;   // update phiByVoidfraction, at least excecute once.
                        }

                    } // end non-orthogonal corrector loop

                    phi = voidfractionf*phiByVoidfraction;
                    #include "continuityErrorPhiPU.H"

                    if (modelType=="B" || modelType=="Bfull")
                        U -= rUA*fvc::grad(p) - f/rho*rUA;
                    else
                        U -= voidfraction*rUA*fvc::grad(p) - f/rho*rUA;

                    U.correctBoundaryConditions();
                    fvOptions.correct(U);

                } // end piso loop
            }

            laminarTransport.correct();
            turbulence->correct();
        }// end solveFlow
        else
        {
            Info << "skipping flow solution." << endl;
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

        particleCloud.clockM().stop("Flow");
        particleCloud.clockM().stop("Global");
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
