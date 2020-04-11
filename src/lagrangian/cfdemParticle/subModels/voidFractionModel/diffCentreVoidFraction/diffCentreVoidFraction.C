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

Author 
    Linhan Ge, University of Newcastle
    linhan.ge@gmail.com
\*---------------------------------------------------------------------------*/

#include "error.H"

#include "diffCentreVoidFraction.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(diffCentreVoidFraction, 0);

addToRunTimeSelectionTable
(
    voidFractionModel,
    diffCentreVoidFraction,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
diffCentreVoidFraction::diffCentreVoidFraction
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    voidFractionModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props"))
{
    maxCellsPerParticle_=readLabel(propsDict_.lookup("maxCellsPerParticle"));
    checkWeightNporosity(propsDict_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

diffCentreVoidFraction::~diffCentreVoidFraction()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void diffCentreVoidFraction::setvoidFraction(double** const& mask,double**& voidfractions,double**& particleWeights,double**& particleVolumes,double**& particleV) const
{
    reAllocArrays();

    scalar radius(-1);
    scalar volume(0);
    scalar cellVol(0);
    scalar scaleVol= weight();
    scalar scaleRadius = pow(porosity(),1./3.);

    for(int index=0; index< particleCloud_.numberOfParticles(); index++)
    {
        //if(mask[index][0])
        //{
            //reset
            for(int subcell=0;subcell<cellsPerParticle_[index][0];subcell++)
            {
                particleWeights[index][subcell]=0;
                particleVolumes[index][subcell]=0;
            }
            cellsPerParticle_[index][0]=1.0;
            particleV[index][0]=0;

            //collecting data
            label particleCenterCellID=particleCloud_.cellIDs()[index][0];
            radius = particleCloud_.radius(index);
            volume = 4.188790205*radius*radius*radius*scaleVol;
            radius *= scaleRadius;
            vector positionCenter=particleCloud_.position(index);

            if (particleCenterCellID >= 0)
            {
                labelHashSet hashSett;    // from OpenFOAM

                //determining label and degree of coveredness of cells covered by the particle
                buildLabelHashSet(radius, positionCenter, particleCenterCellID, hashSett);
                  //Info << "completeSize=" << hashSett.size() << ", completeList =\n" << endl;
                  //for(label i=0;i<hashSett.size();i++) Info << " ," << hashSett.toc()[i] << endl;


                //generating list with cell and subcells
                scalar hashSetLength = hashSett.size();
                if (hashSetLength > maxCellsPerParticle_)
                {
                    FatalError<< "big particle algo found more cells ("<< hashSetLength
                              <<") than storage is prepared ("<<maxCellsPerParticle_<<")" << abort(FatalError);
                }
                else if (hashSetLength > 0)
                {
                    cellsPerParticle_[index][0]=hashSetLength;

                    //making sure that the cell containing the center is the first subcell
                    particleCloud_.cellIDs()[index][0]=particleCenterCellID;
                    //deleting the cell containing the center of the particle
                    hashSett.erase(particleCenterCellID);

                    //==========================//
                    //setting the voidfractions

                    // volumefraction of centre use particle centre method
                    particlefractionNext_[particleCenterCellID] += volume/particleCloud_.mesh().V()[particleCenterCellID];
                    particleWeights[index][0] = 1;
                    particleVolumes[index][0] = volume;
                    particleV[index][0] = volume;

                    for(label i=0;i<hashSetLength-1;i++)
                    {
                        label cellI=hashSett.toc()[i];
                        particleCloud_.cellIDs()[index][i+1]=cellI; //adding subcell represenation
                    }

                    // debug
                    if(index==0)
                    {
                        Info << "particle 0 is represented by " << hashSetLength << "cells" << endl;
                    }
                }//end cells found on this proc
            }// end found cells
        //}// end if masked
    }// end loop all particles
    particlefractionNext_.correctBoundaryConditions();

    // bring voidfraction from Eulerian Field to particle array
    for(int index=0; index< particleCloud_.numberOfParticles(); index++)
    {
        // label cellID = particleCloud_.cellIDs()[index][0];

        voidfractions[index][0] = -1;
    }
}

void diffCentreVoidFraction::buildLabelHashSet
(
    const scalar radius,
    const vector position,
    const label cellID,
    labelHashSet& hashSett
)const
{
    hashSett.insert(cellID);
    //Info<<"cell inserted"<<cellID<<endl;
    const labelList& nc = particleCloud_.mesh().cellCells()[cellID];
    forAll(nc,i){
        label neighbor=nc[i];
        if(!hashSett.found(neighbor) && mag(position-particleCloud_.mesh().C()[neighbor])<radius){
            buildLabelHashSet(radius,position,neighbor,hashSett);
        }
    }

}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace Foam
// ************************************************************************* //
