#!/bin/bash

#===================================================================#
# allrun script for testcase as part of test routine 
# run settlingTest CFD part
# Christoph Goniva - Feb. 2011
#===================================================================#

#- source CFDEM env vars
. ~/.bashrc

#- include functions
source $CFDEM_SRC_DIR/lagrangian/cfdemParticle/etc/functions.sh

#--------------------------------------------------------------------------------#
#- define variables
casePath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"
logpath=$casePath
headerText="run_parallel_CFD"
logfileName="log_$headerText"
solverName="kva_interDyMFoam"
nrProcs="8"
machineFileName="none"   # yourMachinefileName | none
debugMode="off"   
#--------------------------------------------------------------------------------#

#- call function to run a parallel CFD-DEM case
parCFDrun $logpath $logfileName $casePath $headerText $solverName $nrProcs $machineFileName $debugMode