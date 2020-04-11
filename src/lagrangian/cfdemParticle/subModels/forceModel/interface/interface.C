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

#include "interface.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(interface, 0);

addToRunTimeSelectionTable
(
    forceModel,
    interface,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
interface::interface
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    forceModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    alphaFieldName_(propsDict_.lookup("alphaFieldName")),
    alpha_(sm.mesh().lookupObject<volScalarField> (alphaFieldName_)),
    sigmaKFieldName_(propsDict_.lookup("sigmaKFieldName")),
    sigmaK_(sm.mesh().lookupObject<volScalarField> (sigmaKFieldName_)),
    theta_(readScalar(propsDict_.lookup("theta"))),
    alphaLower_(readScalar(propsDict_.lookup("alphaLower"))),
    alphaUpper_(readScalar(propsDict_.lookup("alphaUpper"))),
    backwardInterpolation_(false),
    limitForce_(false),
    C_(readScalar(propsDict_.lookup("C"))),   
    interpolation_(false),
    verbose(false)
{
    // init force sub model
    setForceSubModels(propsDict_);

    // define switches which can be read from dict
    forceSubM(0).setSwitchesList(0,true); // activate treatExplicit switch
    forceSubM(0).setSwitchesList(4,true); // activate search for interpolate switch

    forceSubM(0).setSwitches(0,true);
    // read those switches defined above, if provided in dict
    forceSubM(0).readSwitches();
    //for (int iFSub=0;iFSub<nrForceSubModels();iFSub++)
    //    forceSubM(iFSub).readSwitches();

    particleCloud_.checkCG(false);

    if(propsDict_.found("backwardInterpolation"))  
        backwardInterpolation_ = true;

    if(propsDict_.found("limitForce"))  
        limitForce_ = true;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

interface::~interface()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void interface::setForce() const
{
    volVectorField gradAlpha_ = fvc::grad(alpha_);

    #include "resetAlphaInterpolator.H"
    #include "resetGradAlphaInterpolator.H"
    #include "resetSigmaKInterpolator.H"

    for(int index = 0;index <  particleCloud_.numberOfParticles(); ++index)
    {
        //if(mask[index][0])
        //{
            // definition of spherical particle
            scalar dp = 2*particleCloud_.radius(index);
            vector position = particleCloud_.position(index);
            label cellI = particleCloud_.cellIDs()[index][0];

            if(cellI >-1.0) // particle found on proc domain
            {
                scalar alphap(0);
                vector gradAlpha(0,0,0);
                scalar sigmaK(0);

                if(forceSubM(0).interpolation()) // use intepolated values for alpha (normally off!!!)
                {
                    // make interpolation object for alpha
                    alphap = alphaInterpolator_().interpolate(position,cellI);
                    // make interpolation object for grad(alpha)/|grad(alpha)|
                    gradAlpha = gradAlphaInterpolator_().interpolate(position,cellI);
                    sigmaK = sigmaKInterpolator_().interpolate(position,cellI);
                }
                else if (backwardInterpolation_)
                {
                    vector totalGradAlphaVol(0,0,0);
                    scalar totalSigmaKVol(0);
                    scalar tolVol(0);
 
                    for(int subCell=0;subCell<particleCloud_.cellsPerParticle()[index][0];subCell++) 
                    {
                        label subCellID = particleCloud_.cellIDs()[index][subCell];
                        totalGradAlphaVol += gradAlpha_[subCellID]*particleCloud_.mesh().V()[subCellID];
                        totalSigmaKVol += sigmaK_[subCellID]*particleCloud_.mesh().V()[subCellID];
                        tolVol += particleCloud_.mesh().V()[subCellID];
                    }
                    gradAlpha = totalGradAlphaVol /tolVol;
                    sigmaK = totalSigmaKVol /tolVol;
                }
                else // use cell centered values for alpha
                {
                    //// for any reason fvc::grad(alpha_) cannot be executed here!?
                    //volVectorField gradAlpha=fvc::grad(alpha_);
                    //volVectorField a = gradAlpha/
                    //                   max(mag(gradAlpha),dimensionedScalar("a",dimensionSet(0,-1,0,0,0), SMALL));
                    //n = a[cellI];

                    alphap = alpha_[cellI];
                    gradAlpha = gradAlpha_[cellI];
                    sigmaK = sigmaK_[cellI];
                }

                // Initialize an interfaceForce vector
                vector interfaceForce = Foam::vector(0,0,0);
                scalar Vs = particleCloud_.particleVolume(index);
                scalar rhop = particleCloud_.Density(index);
                
                // Calculate the interfaceForce (range of alphap needed for stability)

                if ( alphaLower_ < alphap && alphap < alphaUpper_ && !limitForce_)
                {
                    interfaceForce = -Vs*sigmaK*gradAlpha*cos(theta_);
                }

                if ( alphaLower_ < alphap && alphap < alphaUpper_ && limitForce_)
                {
                    interfaceForce = C_*Vs*rhop*gradAlpha;
                }

                if(verbose && mag(interfaceForce) > 0)
                {
                Info << "dp = " << dp << endl;
                Info << "position = " << position << endl;
                Info << "cellI = " << cellI << endl;
                Info << "alpha cell = " << alpha_[cellI] << endl;
                Info << "alphap = " << alphap << endl;
                Info << "gradAlpha = " << gradAlpha << endl;
                Info << "interfaceForce = " << interfaceForce << endl;
                Info << "mag(interfaceForce) = " << mag(interfaceForce) << endl;
                }

                // limit interface force
                /*scalar rhoP=3000;
                scalar mP=dp*dp*dp*3.1415/4*rhoP;
                scalar fMax=5*mP*9.81;
                if(mag(interfaceForce)>fMax){
                    interfaceForce /= mag(interfaceForce)/fMax;
                    Info << "interface force is limited to " << interfaceForce << endl;
                }*/

               // write particle based data to global array
               forceSubM(0).partToArray(index,interfaceForce,vector::zero);

               forceSubM(0).passInterfaceForce(index,interfaceForce);

            } // end if particle found on proc domain
        //}// end if in mask
    }// end loop particles
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
