#!/bin/bash

#===================================================================#
# allrun script for testcase as part of test routine 
# run ErgunTestMPI
# Christoph Goniva - Sept. 2010
#===================================================================#

#- define variables
casePath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"

# cd $casePath
# ./parDEMrun.sh

cp $casePath/CFD/0/alpha.water.org $casePath/CFD/0/alpha.water

cd $casePath/CFD
setFields

cd $casePath/CFD/system
cp controlDict_CFDEM controlDict

cd $casePath
./parCFDDEMrun.sh

cd $casePath/CFD
reconstructPar
