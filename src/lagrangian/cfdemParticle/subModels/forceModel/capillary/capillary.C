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

#include "capillary.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(capillary, 0);

addToRunTimeSelectionTable
(
    forceModel,
    capillary,
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
capillary::capillary
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    forceModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    primaryPhaseFieldName_(propsDict_.lookup("primaryPhaseFieldName")),
    alpha_(sm.mesh().lookupObject<volScalarField> (primaryPhaseFieldName_)),
    // velFieldName_(propsDict_.lookup("velFieldName")),
    // U_(sm.mesh().lookupObject<volVectorField> (velFieldName_)),
    sigma_(readScalar(propsDict_.lookup("sigma"))),
    theta_(readScalar(propsDict_.lookup("theta"))),
    alphaThreshold_(readScalar(propsDict_.lookup("alphaThreshold"))),
    deltaAlphaIn_(readScalar(propsDict_.lookup("deltaAlphaIn"))),
    deltaAlphaOut_(readScalar(propsDict_.lookup("deltaAlphaOut"))),
    C1_(1.0),
    C2_(1.0),
    backwardInterpolation_(false),
    interpolation_(false),
    verbose(false)
    // centers(sm.mesh().C()),
    // vols(sm.mesh().V())
    /*alpha05Dict
    (
        IOobject
        (
            "alpha05Dict",
            sm.mesh().time().constant(),
            sm.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    alpha05
    (
        "alpha05",
        sm.mesh(),
        alpha05Dict
    )*/
{
    if (propsDict_.found("C1")) C1_=readScalar(propsDict_.lookup("C1"));
    if (propsDict_.found("C2")) C2_=readScalar(propsDict_.lookup("C2"));
    
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
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

capillary::~capillary()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void capillary::setForce() const
{
    // alpha05.update();

    volVectorField gradAlpha_ = fvc::grad(alpha_);
    #include "resetAlphaInterpolator.H"
    #include "resetGradAlphaInterpolator.H"

    // vector Us(0,0,0);

    // const pointField&  alpha05point = alpha05.points();   //points on the surface
    // faceList surfaceFace = alpha05.faces();               //faces on the surface
    // pointField centre = alpha05point.Cf();

    // List<vector> centre(surfaceFace.size(),vector(0,0,0));
    
    // pointField pp(alpha05point.size(),vector(0,0,0));

    /*forAll(surfaceFace,facei)
    {    
        // const face& ll = surfaceFace[facei];  //list of  point label of facei
        // const pointField& pll = ll.points(alpha05point);
        // centre[facei] = ll.centre(pll); 
        // Info << centre[facei] << endl;
        Info << alpha05.Cf()[facei] << endl;
    }*/
    //forAll(alpha05point,iP)
    //{ 
        // points[iP] = vector(_coo[iP].x, _coo[iP].y, _coo[iP].z);
        // Sout << alpha05point[iP] << endl;
        
    //}

    /*label pointNo = alpha05point.size();
    reduce(pointNo,sumOp<label>());
    Info << pointNo << endl;*/

    // sum of U by volume
    /*scalar sumUX = 0;
    scalar sumUY = 0;
    scalar sumUZ = 0;
    // sum of volume
    scalar sumVol = 0;

    forAll(centers, cellI)
    {
        if ( alpha_[cellI] < alphaThreshold_)
        {
            //Rising Velocity
            sumUX   += (1 - alpha_[cellI])*vols[cellI]*U_[cellI].x();
            sumUY   += (1 - alpha_[cellI])*vols[cellI]*U_[cellI].y();
            sumUZ   += (1 - alpha_[cellI])*vols[cellI]*U_[cellI].z();
            sumVol  += (1 - alpha_[cellI])*vols[cellI];
        }
    }
    
    reduce(sumUX, sumOp<scalar>());
    reduce(sumUY, sumOp<scalar>());
    reduce(sumUZ, sumOp<scalar>());
    reduce(sumVol, sumOp<scalar>());

    // rising velocity
    scalar vx = sumUX/max(sumVol, SMALL);
    scalar vy = sumUY/max(sumVol, SMALL);
    scalar vz = sumUZ/max(sumVol, SMALL);*/

    for(int index = 0;index <  particleCloud_.numberOfParticles(); ++index)
    {
        //if(mask[index][0])
        //{
            // definition of sphericalparticle
            scalar dp = 2*particleCloud_.radius(index);
            vector position = particleCloud_.position(index);
            
            /*List<vector> pointOnLine;
            int nuOfPoints = 1e3;
            pointOnLine.setSize(nuOfPoints, vector(0,0,0)));
            
            for (int i=1, i <= nuOfPoints, i++)
            {
                pointOnLine(i) = ((cg - position)/1e6*i + position);    
            }*/ 

            label cellI = particleCloud_.cellIDs()[index][0];
            // vector interfacePoint(0,0,0);
            // scalar h(0);

            if(cellI >-1.0) // particle found on proc domain
            {
                scalar alphap;
                vector n;
                vector gradAlphap(0,0,0);
                // Us = particleCloud_.velocity(index);
                /*for (int i = 1, i<=nuOfPoints,i++)
                {
                    alphaPoint = alphaInterpolator_().interpolate(pointOnLine(i),cellI);
                    if abs( alphaPoint - alphaThreshold_ ) < 1e6
                    {
                        interfacePoint = pointOnLine(i);
                    }
                }
                // separate distance between the point and the interface
                h = sqrt((position.x()-interfacePoint.x())*(position.x()-interfacePoint.x())
                   +(position.y()-interfacePoint.y())*(position.y()-interfacePoint.y())
                   +(position.z()-interfacePoint.z())*(position.z()-interfacePoint.z()));*/

                if(forceSubM(0).interpolation()) // use intepolated values for alpha (normally off!!!)
                {
                    // make interpolation object for alpha
                    alphap = alphaInterpolator_().interpolate(position,cellI);

                    // make interpolation object for grad(alpha)/|grad(alpha)|
                    gradAlphap = gradAlphaInterpolator_().interpolate(position,cellI);
                    n = gradAlphap/max(mag(gradAlphap),SMALL);
                }
                else if (backwardInterpolation_)
                {
                    vector totalGradAlphaVol(0,0,0);
                    scalar tolVol(0);
 
                    for(int subCell=0;subCell<particleCloud_.cellsPerParticle()[index][0];subCell++) 
                    {
                        label subCellID = particleCloud_.cellIDs()[index][subCell];
                        totalGradAlphaVol += gradAlpha_[subCellID]*particleCloud_.mesh().V()[subCellID];
                        tolVol += particleCloud_.mesh().V()[subCellID];
                    }
                    gradAlphap = totalGradAlphaVol /tolVol;
                    n = gradAlphap/max(mag(gradAlphap),SMALL);
                }
                else // use cell centered values for alpha
                {
                    //// for any reason fvc::grad(alpha_) cannot be executed here!?
                    //volVectorField gradAlpha=fvc::grad(alpha_);
                    //volVectorField a = gradAlpha/
                    //                   max(mag(gradAlpha),dimensionedScalar("a",dimensionSet(0,-1,0,0,0), SMALL));
                    //n = a[cellI];

                    alphap = alpha_[cellI];
                    gradAlphap = gradAlpha_[cellI];
                    n = gradAlphap/max(mag(gradAlphap),SMALL);
                }

                // velocity at gradient direction
                // vector Ub(vx,vy,vz); 

                // bubble rising velocity along the gradient direction                
                // vector Ubn = Ub & n/(n & n) * n;
                // relative velocity at gradient
                // vector Usn = Us & n/(n & n) * n;
                // relative velocity along the gradient direction
                // vector Urn = Ubn-Usn;
                // scalar magUr = mag(Ur);
                // Initialize an capillaryForce vector
                vector capillaryForce = Foam::vector(0,0,0);
                
                // Calculate the capillaryForce (range of alphap needed for stability)
                // Calculate estimate attachment force as
                // |6*sigma*sin(pi-theta/2)*sin(pi+theta/2)|*2*pi*dp
                if ((alphaThreshold_-deltaAlphaIn_) < alphap && alphap < (alphaThreshold_+deltaAlphaOut_))
                {
                scalar Fatt = 0.5*M_PI*dp*sigma_*(1-cos(theta_));
                // C_ can be specified as the maximum value as 1/(dacayFactor*sqrt(2*pi)) to realize the Gaussian distribution
                // be careful with the sign, gradient points from low pressure to high pressure
                // C2 is a damping coefficient
                capillaryForce = -1*n*Fatt*C1_+(alphaThreshold_+deltaAlphaOut_-alphap)/(deltaAlphaIn_+deltaAlphaOut_)*n*Fatt*C2_;
                                // -1*n*Fatt*sin(alphap*M_PI)*C1_+Urn*C2_;
                                // -1*n*Fatt*C_*sin(alphap*M_PI);
                                // * exp(-0.5*((alphap-alphaCentre_)/decayFactor_)*((alphap-alphaCentre_)/decayFactor_));
                }

                if(verbose && mag(capillaryForce) > 0)
                {
                Info << "dp = " << dp << endl;
                Info << "position = " << position << endl;
                Info << "cellI = " << cellI << endl;
                Info << "alpha cell = " << alpha_[cellI] << endl;
                Info << "alphap = " << alphap << endl;
                Info << "n = " << n << endl;
                Info << "capillaryForce = " << capillaryForce << endl;
                Info << "mag(capillaryForce) = " << mag(capillaryForce) << endl;
                }

                // limit capillary force
                /*scalar rhoP=3000;
                scalar mP=dp*dp*dp*3.1415/4*rhoP;
                scalar fMax=5*mP*9.81;
                if(mag(capillaryForce)>fMax){
                    capillaryForce /= mag(capillaryForce)/fMax;
                    Info << "capillary force is limited to " << capillaryForce << endl;
                }*/

               // write particle based data to global array
               forceSubM(0).partToArray(index,capillaryForce,vector::zero);

            } // end if particle found on proc domain
        //}// end if in mask
    }// end loop particles
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// calculate the centre of the bubble
    /*const volVectorField& centers = mesh().C();
    scalar sumX = 0;
    scalar sumY = 0;
    scalar sumZ = 0;
    
    const objectRegistry& db = mesh().thisDb();
    const dictionary& transportProperties = 
          db.lookupObject<IOdictionary>
          (
            "transportProperties"
          );

    dictionary alpha2
    (
        transportProperties.subDict("air")
    );
        
    dimensionedScalar rhoBubble
    (
        "rhoBubble",
        dimDensity,
        alpha2.lookup("rho")
    );

    forAll(centers, I)
    {
        if ( alpha_[I] < th && centers[I].y() < ht )
        {
            //Mass center
            sumY    += (1 - alpha_[I])*rhoBubble.value()*vols[I]*centers[I].y();
            sumX    += (1 - alpha_[I])*rhoBubble.value()*vols[I]*centers[I].x();
            sumZ    += (1 - alpha_[I])*rhoBubble.value()*vols[I]*centers[I].z();
            sumMass += (1 - alpha_[I])*rhoBubble.value()*vols[I];
        }
    }
    
    reduce(sumY, sumOp<scalar>());
    reduce(sumX, sumOp<scalar>());
    reduce(sumZ, sumOp<scalar>());
        
    scalar cgX = sumX/max(sumMass, SMALL);
    scalar cgY = sumY/max(sumMass, SMALL);
    scalar cgZ = sumZ/max(sumMass, SMALL);

    vectro cg(cgX,cgY,cgZ);*/
// ************************************************************************* //
