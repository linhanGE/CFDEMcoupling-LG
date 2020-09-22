#!/bin/bash

#===================================================================#
# allrun script for testcase as part of test routine 
# run ErgunTestMPI
# Christoph Goniva - Sept. 2010
#===================================================================#

#- define variables
casePath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"

# $casePath/parDEMrun.sh

cd $casePath

./parCFDDEMrun.sh

cd $casePath/CFD

reconstructPar

cd $casePath/CFD

rm -rf process*

# foamToVTK -fields '(p U alpha.water)' -excludePatches '(walls p1 p2)'

# rm -rf $casePath/CFD/0.*

# rm -rf $casePath/DEM/post/*runboundingBox*.vtk
