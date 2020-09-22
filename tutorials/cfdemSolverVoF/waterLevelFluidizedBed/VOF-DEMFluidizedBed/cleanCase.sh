#!/bin/bash

#===================================================================#
# allrun script for testcase as part of test routine 
# run settlingTest CFD part
# Christoph Goniva - Feb. 2011
#===================================================================#

#- include functions
casePath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"

rm -rf $casePath/CFD/0.*
rm -rf $casePath/CFD/1*
rm -rf $casePath/CFD/2*
rm -rf $casePath/CFD/3*
rm -rf $casePath/CFD/4*
rm -rf $casePath/CFD/5*
rm -rf $casePath/CFD/6*
rm -rf $casePath/CFD/7*
rm -rf $casePath/CFD/8*
rm -rf $casePath/CFD/9*

rm -rf $casePath/CFD/log*

rm -rf $casePath/CFD/post*

rm -rf $casePath/CFD/proce*

rm -rf $casePath/DEM/post/dump*

rm -rf $casePath/DEM/post/*.txt

rm -rf $casePath/DEM/post/restart/*.*

rm -rf $casePath/DEM/post/*.vtk

rm -rf $casePath/DEM/post/*.run

rm -rf $casePath/DEM/log.*

rm -rf $casePath/DEM/post/gp/*.local

rm -rf $casePath/DEM/post/gp/*.txt 

rm -rf $casePath/DEM/post/mfp/*.txt

rm -rf $casePath/DEM/post/v/*.txt

rm -rf $casePath/log*

rm -rf $casePath/*.log

rm -rf $casePath/bubbleData

rm -rf $casePath/*.sh.*

rm -rf $casePath/CFD/0/alpha.water

rm -rf $casePath/log*

rm -rf $casePath/bubbleData

rm -rf $casePath/CFD/dynamicCode