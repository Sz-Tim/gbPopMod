### Szewczyk TM, T Lee, MJ Ducey, ME Aiello-Lammens, H Bibaud, JM Allen. 2019. Local management in a regional context: Simulations with process-based species distribution models. *Ecological Modelling* 413(1): 108827  

---  

Code repository accompanying [Szewczyk et al. 2019](https://doi.org/10.1016/j.ecolmodel.2019.108827)  

---  

This R package runs a population-level demographic cellular automata model. The life history structure and default parameterization are based on [glossy buckthorn](https://en.wikipedia.org/wiki/Frangula_alnus) (*Frangula alnus*) in New England. The model is spatially explicit and temporally dynamic, simulating the growth and spread of a species on a gridded landscape. Population and life history dynamics occur within each occupied grid cell, driven by vital rates and demographic parameters which can be global or dependent on the cell environment. Cells are connected through mechanistic short distance dispersal (SDD) and rare long distance dispersal (LDD). Management actions can target specific life stages, such as reducing seedling establishment (e.g., through planted ground cover) or reducing the number of adults (e.g., through cutting and/or spraying).  In each year, the simulation occurs according to the following series of steps:  

![Simulation structure. In each time step, the sequence begins with fruit production based on the number of adults in a given cell. Boxes represent containers, with associated abundances, and text along arrows represent parameters. Parameters shown in blue vary by land cover type, while parameters in black are global.](sensitivity/model_outline.jpeg)

with the following parameters and containers for each cell *i*, each time *t*, and juvenile ages *k* = `1:(m-1)`:  

| Symbol | Quantity |
| :--- | :--- |
| *N<sub>it</sub>* | Number of adults |
| *M<sub>itk</sub>* | Number of juveniles |
| *F<sub>it</sub>* | Number of fruits |
| *S<sub>it</sub>* | Number of seeds deposited |
| *B<sub>it</sub>* | Number of seeds in seed bank |
| *f<sub>i</sub>* | Mean flowering probability |
| *μ<sub>i</sub>* | Mean number of fruits from flowering adult |
| *γ* | Mean number of seeds per fruit |
| *m<sub>i</sub>* | Age at adulthood |
| *c<sub>i</sub>* | Pr(fruit is consumed by a bird) |
| *r* | Rat parameter for SDD exponential kernel (units: cells) |
| *sdd<sub>max</sub>* | Maximum SDD distance (units: cells) |
| *J<sub>i</sub>* | Number of SDD neighbors (calculated from *sdd<sub>max</sub>*) |
| *δ<sub>ji</sub>* | Pr(emigrant from *j* is deposited in *i*) |
| *n<sub>ldd</sub>* | Number of LDD events |
| *s<sub>c</sub>* | Pr(seed survival \| consumed) |
| *s<sub>B</sub>* | Pr(seed survival in the seed bank) |
| *s<sub>M,i</sub>* | Pr(juvenile survival) |
| *s<sub>N,i</sub>* | Pr(adult survival) |
| *K<sub>i</sub>* | Carrying capacity for adults |
| *g<sub>D</sub>* | Pr(germinating directly: same year as produced) |
| *g<sub>B</sub>* | Pr(germinating from the seed bank) |
| *p<sub>i</sub>* | Pr(establishing once germinated) |

Note that in the code, the ages `1:(m-1)` are included as layers in the object `N` rather than as a separate object `M`.

Parameters that vary with among cells can be estimated through regressions with environmental covariates [(Merow et al 2017)](https://doi.org/10.1073/pnas.1609633114) as appropriate for the species and system. Density dependence can be implemented either as a hard cap on the abundance within each cell [(Merow et al 2011)](https://doi.org/10.1086/660295) or through parameters such as the establishment probability, survival probabilities, or fruiting probabilities [(Ellner & Rees 2006)](https://doi.org/10.1086/499438).

# Running the model  
The model is run using the function `run_sim()`, though some set up is necessary for initializing the landscape. The script `vignettes/management_controls.R` assigns parameter values, establishes management strategies, and builds the landscape before running a set number of stochastic simulations. The script `hpc/hpc_wrap.R` runs a global sensitivity analysis, where many parameters are varied simultaneously. 
