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
    cfdemSolverVoF

Description
    Transient solver for incompressible gas-liquid-solid flow.
    Turbulence modelling is generic, i.e. laminar, RAS or LES may be selected.
    The code is an evolution of the solver pimpleFoam in OpenFOAM(R) 1.6,
    where additional functionality for CFD-DEM coupling is added.
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "CMULES.H"
#include "EulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"

#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "turbulentTransportModel.H"
#include "pisoControl.H"
#include "fvOptions.H"
#include "CorrectPhi.H"
#include "fixedFluxPressureFvPatchScalarField.H"
#include "fvcSmooth.H"

#include "cfdemCloud.H"
#include "explicitCouple.H"
#include "clockModel.H"
#include "smoothingModel.H"
#include "forceModel.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"
	#include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createTimeControls.H"
	#include "initContinuityErrs.H"
    #include "createFields.H"
	#include "createAlphaFluxes.H"
    #include "createFvOptions.H"
    #include "correctPhi.H"
	
	turbulence->validate();

	// create cfdemCloud
	
	cfdemCloud particleCloud(mesh);

	#include "checkModelType.H"
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info<< "\nStarting time loop\n" << endl;
	
    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;
		
        #include "readTimeControls.H"
		#include "CourantNo.H"
        #include "alphaCourantNo.H"
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

		voidfractionf = fvc::interpolate(voidfraction);
		voidfractionPhic = voidfractionf*phi;
		
        //Force Checks
        #include "forceCheckEx.H" 

        particleCloud.clockM().stop("Coupling");

        particleCloud.clockM().start(26,"Flow");
		
		if(particleCloud.solveFlow())
		{
			#include "alphaControls.H"
			#include "alphaEqnSubCycle.H"
		 
			mixture.correct();
			voidfractionRho = voidfraction*rho;

			#include "UEqn.H"

			while (piso.correct())
			{
				#include "pEqn.H"
			} 
		
			turbulence->correct();
		}
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
