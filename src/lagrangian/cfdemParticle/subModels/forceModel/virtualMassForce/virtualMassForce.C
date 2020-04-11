/*---------------------------------------------------------------------------*\
    CFDEMcoupling - Open Source CFD-DEM coupling

    CFDEMcoupling is part of the CFDEMproject
    www.cfdem.com
                                Christoph Goniva, christoph.goniva@cfdem.com
                                Copyright 2009-2012 JKU Linz
                                Copyright 2012-     DCS Computing GmbH, Linz
-------------------------------------------------------------------------------
License
    This file is part of CFDEMcoupling.

    CFDEMcoupling is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    CFDEMcoupling is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with CFDEMcoupling; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description
    This code is designed to realize coupled CFD-DEM simulations using LIGGGHTS
    and OpenFOAM(R). Note: this code is not part of OpenFOAM(R) (see DISCLAIMER).
\*---------------------------------------------------------------------------*/

#include "error.H"

#include "virtualMassForce.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

#define NOTONCPU 9999

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(virtualMassForce, 0);

addToRunTimeSelectionTable
(
    forceModel,
    virtualMassForce,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
virtualMassForce::virtualMassForce
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    forceModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    velFieldName_(propsDict_.lookup("velFieldName")),
    U_(sm.mesh().lookupObject<volVectorField> (velFieldName_)),
    Cadd_(0.5),
    backwardInterpolation_(false)
{
    // init force sub model
    setForceSubModels(propsDict_);
    // define switches which can be read from dict
    forceSubM(0).setSwitchesList(0,true); // activate treatExplicit switch
    forceSubM(0).setSwitchesList(1,true); // activate treatForceDEM switch (DEM side only treatment)
    forceSubM(0).setSwitchesList(4,true); // activate search for interpolate switch
    

    for (int iFSub=0;iFSub<nrForceSubModels();iFSub++)
        forceSubM(iFSub).readSwitches();

    Info << "WARNING: be sure that LIGGGHTS integration takes ddtv_p implicitly into account! \n";

    if(propsDict_.found("Cadd"))
    {
        Cadd_ = readScalar(propsDict_.lookup("Cadd")); 
        Info << "Virtual mass model: using non-standard Cadd = " << Cadd_ << endl;
    }

    particleCloud_.checkCG(true);

    //Append the field names to be probed
    particleCloud_.probeM().initialize(typeName, typeName+".logDat");
    particleCloud_.probeM().vectorFields_.append("virtualMassForce"); //first entry must the be the force
    particleCloud_.probeM().vectorFields_.append("DDtU");
    particleCloud_.probeM().scalarFields_.append("Vs");
    particleCloud_.probeM().scalarFields_.append("rho");
    particleCloud_.probeM().writeHeader();

    if (propsDict_.found("backwardInterpolation")) 
    {
        backwardInterpolation_=true;
        Info << "Use average backward interpolation" <<endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

virtualMassForce::~virtualMassForce()
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void virtualMassForce::setForce() const
{

    vector position(0,0,0);
    vector Ufluid(0,0,0);
    vector DDtU(0,0,0);
    
    //Compute extra vfields in case it is needed
    volVectorField DDtU_(0.0*U_/U_.mesh().time().deltaT());
    surfaceScalarField phi_ = fvc::flux(U_);
    DDtU_ = fvc::DDt(phi_, U_); //Total Derivative of fluid velocity

    #include "resetUInterpolator.H"
    #include "resetDDtUInterpolator.H"
    #include "setupProbeModel.H"

    for(int index = 0;index <  particleCloud_.numberOfParticles(); index++)
    {
            vector virtualMassForce(0,0,0);
            label cellI = particleCloud_.cellIDs()[index][0];

            if (cellI > -1) // particle Found
            {

                if(forceSubM(0).interpolation()) 
                {
	                position     = particleCloud_.position(index);
                    Ufluid       = UInterpolator_().interpolate(position,cellI);
                    DDtU = DDtUInterpolator_().interpolate(position,cellI);
                }
                else if (backwardInterpolation_)
                {
                    vector totalUfluidVol(0,0,0);
                    vector totalDDtUVol(0,0,0);
                    scalar tolVol(0);

                    for(int subCell=0;subCell<particleCloud_.cellsPerParticle()[index][0];subCell++) 
                    {
                        label subCellID = particleCloud_.cellIDs()[index][subCell];
                        totalUfluidVol += U_[subCellID]*particleCloud_.mesh().V()[subCellID];
                        totalDDtUVol += DDtU_[subCellID]*particleCloud_.mesh().V()[subCellID];
                        tolVol += particleCloud_.mesh().V()[subCellID];
                    }
                    Ufluid = totalUfluidVol /tolVol;
                    DDtU = totalDDtUVol /tolVol;
                }
                else
                {
                    Ufluid = U_[cellI];
                    DDtU = DDtU_[cellI];
                }
                
                scalar ds = 2*particleCloud_.radius(index);
                scalar Vs = ds*ds*ds*M_PI/6;
                scalar rho = forceSubM(0).rhoField()[cellI];
                
                virtualMassForce = Cadd_ * rho * Vs * DDtU;

                if( forceSubM(0).verbose() ) //&& index>100 && index < 105)
                {
                    Pout << "index / cellI = " << index << tab << cellI << endl;
                    Pout << "position = " << particleCloud_.position(index) << endl;
                }

                //Set value fields and write the probe
                if(probeIt_)
                {
                    #include "setupProbeModelfields.H"
                    // Note: for other than ext one could use vValues.append(x)
                    // instead of setSize
                    vValues.setSize(vValues.size()+1, virtualMassForce);           //first entry must the be the force
                    vValues.setSize(vValues.size()+1, DDtU);
                    sValues.setSize(sValues.size()+1, Vs);
                    sValues.setSize(sValues.size()+1, rho);
                    particleCloud_.probeM().writeProbe(index, sValues, vValues);
                }
            }
    
            // write particle based data to global array
            forceSubM(0).partToArray(index,virtualMassForce,vector::zero);
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
