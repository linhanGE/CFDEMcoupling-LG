# This is a modified version of [CFDEM](https://www.cfdem.com/).
# Author: Linhan Ge, School of Engineering, the University of Newcastle.
# Typical features:
## CFD-DEM solver with diffusion averaging method developed by [Sun and Xiao](https://www.sciencedirect.com/science/article/pii/S030193221500186X)
* Grid cells that smaller than particles can be used to resolve the fluid flow.
* Details can be found in the below paper:

[Ge, L., Evans, G. M., & Moreno-Atanasio, R. (2020). CFD-DEM investigation of the interaction between a particle swarm and a stationary bubble: Particle-bubble collision efficiency. Powder Technology.](https://www.sciencedirect.com/science/article/pii/S0032591020302102)


 *Please cite my paper if you use the solver in your research.*

## VOF-DEM solver with diffusion averaging method developed by [Sun and Xiao](https://www.sciencedirect.com/science/article/pii/S030193221500186X)
* The solver is a gas-liquid-solid three phase flow solver which is capable of resolving gas-liquid interfaces using extremely fine meshes whilst guaranteeing physical particle dynamics.
* Details can be found in the below paper:


[Ge, L., Peng, Z., Moreno-Atanasio, R., Doroodchi, E., & Evans, G. M. (2020). A Three-dimensional VOF-DEM Model for Simulating Particle Dynamics in the Liquid Slugs of a Vertical Gas-Liquid-Solid Taylor Flow Microreactor. Industrial & Engineering Chemistry Research.](https://pubs.acs.org/doi/abs/10.1021/acs.iecr.0c00108)


*Please cite my paper if you use the solver in your research.*

# Installation procedures:
* Install [OpenFOAM 5.x](https://openfoamwiki.net/index.php/Installation/Linux/OpenFOAM-5.x/Ubuntu).
* Install the smoothing curvature feature for interFoam developed by [Kevin van As](https://github.com/floquation/OF-kva_interfaceProperties).
* Install CFDEM as described in the [CFDEM manual](https://www.cfdem.com/media/CFDEM/docu/CFDEMcoupling_Manual.html) but using the CFDEM repository here. Noting that the folder name should be changed from *CFDEMcoupling-LG* to *CFDEMcoupling-PUBLIC-5.x*.
* Please set up your own cases following the examples in the tutorials.
# List of new solvers
* cfdemSolverDiffusion (CFD-DEM)
* cfdemSolverVoF (VOF-DEM)
* cfdemSolverVoFCyclic (VOF-DEM compatible with periodic domains)
