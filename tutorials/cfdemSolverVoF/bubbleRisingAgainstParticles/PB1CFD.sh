#!/bin/bash
# 
# Script to run CFDEM on Linux Grid (v1)
#
########################################################
#PBS -l select=1:ncpus=8:mem=8gb
#PBS -l walltime=80:00:00
#PBS -k oe
#PBS -M c3216945@uon.edu.au
#PBS -m bae

# Edit where my simulation is and the name of the file
SIMDIR="/home/c3216945/paperPB/PB1"

source /home/c3216945/OpenFOAM/OpenFOAM-5.x/etc/bashrc
source /home/c3216945/CFDEM/CFDEMcoupling-PUBLIC-5.x/src/lagrangian/cfdemParticle/etc/bashrc

cd $SIMDIR

./runCFD.sh

exit 0


