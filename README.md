# Gravity driven infiltration with the development of a saturated region
## Authors
- Mohammad Afzal Shadab (mashadab@utexas.edu)
- Marc Andre Hesse (mhesse@jsg.utexas.edu)

## Affiliation
Oden Institute for Computational Engineering and Sciences / Jackson School of Geosciences / University of Texas Institute for Geophysics
The University of Texas at Austin

## 5-line Summary
Hyperbolic (kinematic wave) analysis in infiltration theory is a technique to estimate the ponding time and runoff estimates in low-textured soils. The proposed theory can be applied to general soil profiles. The codes provide the analytical solutions for infiltration in soils with two-layer, exponential and power-law porosity decay profiles. The results have been validated with numerical solver and Hydrus simulations. The codes given correspond to the figures from the research paper. 

## Citation
[1] Shadab, M.A. and Hesse, M.A., 202X. Gravity driven infiltration with the development of a saturated region, Water Resources Research (in preparation).
[2] Shadab, M.A. and Hesse, M.A., 2021, December. Fluid Infiltration in Unsaturated Porous Medium with The Development of a Saturated Region. In AGU Fall Meeting 2021. AGU.

## Getting Started

### Overview

Understanding the controls on the infiltration of precipitation into soil is an important problem in hydrology and its representation in heterogeneous media is still challenging. Here we investigate the reduction of infiltration by the development of a saturated region due to the downward decrease in porosity and permeability. The formation of a saturated region leads to the back-filling of the pore space above, so that the saturated region expands both upward and downward. Practically this problem in important because soil porosity commonly declines downward and leads to a reduction of infiltration. We consider gravity-driven infiltration using the hyperbolic limit of Richards' equation. The analytical results for all soil profiles are validated with numerical simulations and moreover two-layered profile is also validated with Hydrus-1D results. This approach helps analyze infiltration in more general low-textured soil profiles at different rainfall rates and provides ponding time and runoff estimates. Lastly, solving the hyperbolic-elliptic problem is also computationally more efficient than solving the full Richards' equation in the limit of negligible capillary forces, because the non-linearity can be integrated explicitly. 

Key points of the present work:
1. Derived relations for infiltration due to transitional rainfall in a dry soil with general porosity decay with depth
2. Integrated the fully-saturated elliptic region inside an otherwise unsaturated hyperbolic region

<p align="center">
<img src="./Figures/Cover.png" height="370">
</p>
Figure : Schematic showing gravity-driven infiltration in a soil with porosity ?? decay with depth. The colors blue, brown, white refer to water, soil and gas respectively. A fully saturated region ??(t) develops within an otherwise unsaturated domain. The saturated-unsaturated region boundary ?????(t) has a boundary condition of water-gas pressure equivalence. The saturated region expands with time as the boundary ?????(t) moves in outward direction.

### Dependences

hyperbolic-infiltration-theory requires the following packages to function:
- [Python](https://www.python.org/) version 3.5+
- [Numpy](http://www.numpy.org/) >= 1.16
- [scipy](https://www.scipy.org/) >=1.5
- [matplotlib](https://matplotlib.org/) >=3.3.4


### Quick Usage
After cloning the repository and installing the required libraries, run the python files corresponding to the figure numbers as given in the paper. Codes can be run either directly or on an IDE such as Anaconda Spyder. `supporting_tools.py` is the file containing the auxillaries. Figure 14 and 15 correspond to a Hydrus 1D file for infiltration in two-layered soil. Use `%matplotlib qt` to change the plotting from inline to separate window in Spyder. The resulting plots are collected in the folder `Figures`.

### Non-dimensionalization
The depth coordinate `z` is scaled with characteristic length `z_0`, time is scaled with characteristic time `z_0/f_c` and infiltration rate `I(t)` (or volumetric flux) is scaled with the infiltration capacity `f_c`. Therefore, the dimensionless variables are `z'=z/z_0`, `t'=tf_c/z_0`, and $`f'=I(t)/f_c`$.
