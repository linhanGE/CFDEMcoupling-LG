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

Class
    capillary

SourceFiles
    capillary.C

\*---------------------------------------------------------------------------*/

#ifndef capillary_H
#define capillary_H

#include "forceModel.H"
#include "sampledIsoSurface.H"
#include "interpolationCellPoint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class capillary Declaration
\*---------------------------------------------------------------------------*/

class capillary
:
    public forceModel
{
private:
    dictionary propsDict_;

    word primaryPhaseFieldName_;

    const volScalarField& alpha_;    // ref to VOF voidfraction

    // word velFieldName_;
    
    // const volVectorField& U_;

    const scalar sigma_;

    const scalar theta_;

    const scalar alphaThreshold_;

    const scalar  deltaAlphaIn_;

    const scalar  deltaAlphaOut_;

    scalar C1_; // model constant

    scalar C2_; // damping coefficient

    bool backwardInterpolation_;

    bool interpolation_; // use interpolated field values (does not work for fine grids???)

    bool verbose;

    // const volVectorField& centers;

    // const scalarField& vols;
    
    // IOdictionary alpha05Dict;

    // mutable sampledIsoSurface alpha05;
public:

    //- Runtime type information
    TypeName("capillary");


    // Constructors

        //- Construct from components
        capillary
        (
            const dictionary& dict,
            cfdemCloud& sm
        );

    // Destructor

        ~capillary();


    // Member Functions
        void setForce() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
