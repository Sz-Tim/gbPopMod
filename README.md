# gbPopMod  

This R package runs a population-level demographic cellular automata model. The life history structure and default parameterization are based on glossy buckthorn (Frangula alnus) in New England. The model is spatially explicit and temporally dynamic, and allows for management targeted to particular life stages in particular grid cells. In each year, the simulation occurs according to the following series of steps:  

![Simulation structure](sensitivity/model_outline.jpeg)

# Running the model  
The model is run using the function `run_sim()`, though some set up is necessary for initializing the landscape. The script `vignettes/management_controls.R` assigns parameter values, establishes management strategies, and builds the landscape before running a set number of stochastic simulations. The script `hpc/hpc_wrap.R` runs a global sensitivity analysis, where many parameters are varied simultaneously. 