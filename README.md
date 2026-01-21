# Overview

Kamal Chandrakar modified the ODT model to include water vapor. The moist RBC convection also accounts for the effect of water vapor on buoyancy in turbulence generation. This model can be used to investigate the mixing process and supersaturation PDF for given boundary conditions  without aerosol or other microphysical processes. The bottom and top boundary conditions determine the water vapor fluxes into the domain and the supersaturation PDF. At a very high resolution, the ODT model can resolve the boundary layer near the bottom and top surfaces. Note that the supersaturation inside the domain will increase and reach a steady state at values larger than those observed since the aerosol injection and microphysical processes are not implemented.  For more information, see Chandarkar et al. (2019)
 
# Contributions
* **Kamal Chandrakar** (v1.0.0) modified the ODT model code to include water vapor as another scalar in addition to temperature. This includes adding required state variables, initial and boundary conditions, outputs for water vapor/supersaturation, and water vapor effects on turbulence via buoyancy.
* **Mani Rajagopal** (v2.0.0) modified the moist-ODT code  from Kamal Chandrakar to output eddy information, scalar probability distribution, covariances, and evolution of domain statistics. Also included is code for input and output to the respective folders, as well as a random initial perturbation seed (as command line argument) to run ensembles of experiments.

# References
Chandrakar, K. K., Cantrell, W., Krueger, S., Shaw, R. A., & Wunsch, S. (2020). Supersaturation fluctuations in moist turbulent Rayleigh–Bénard convection: A two-scalar transport problem. Journal of Fluid Mechanics, 884, A19.
