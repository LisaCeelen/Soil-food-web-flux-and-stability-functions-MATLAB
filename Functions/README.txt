This folder contains the functions created. Functions can be used to compute carbon and nitrogen fluxes and stability of soil food webs of various sizes. 

Data must agree to the following: 
Species must be given in a vertical character list with species names in order of computation. With top predator at first position.
Important to note is that a functional group must be computed after all its predators have been computed, 
In this way the death by predation value is complete. A functional group must always come after all it's predators. 

Biomass, assimilation efficiency, production efficency, natural death rate and C:N ratios are vertical arrays, where indices
agree with the indices defined in the Species list. 

Feeding preferences is a X*X matrix, where X is the number of species in the soil food web, rows represent the predators and columns represent the preys. 
Inactive feeding relationships are 0. Active feeding relationships are only present in top-right half of the matrix. 

