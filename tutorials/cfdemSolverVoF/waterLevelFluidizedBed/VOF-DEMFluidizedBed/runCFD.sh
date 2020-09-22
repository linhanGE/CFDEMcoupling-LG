#!/bin/bash

#===================================================================#
# allrun script for testcase as part of test routine 
# run settlingTest
# Christoph Goniva - Sept. 2010
#===================================================================#

#- define variables
casePath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"

cp $casePath/CFD/0/alpha.water.orig $casePath/CFD/0/alpha.water

cd $casePath/CFD/system

cp -r controlDict_CFD controlDict

cd $casePath/CFD

setFields

. $casePath/CFD/parCFDrun.sh

cd $casePath/CFD

reconstructPar