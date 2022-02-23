# CPEG.jl


A Julia implementation of [CPEG](http://roboticexplorationlab.org/papers/ktracy_cpeg.pdf). An existing implementation of this algorithm is also available in [EntryGuidance.jl](https://github.com/RoboticExplorationLab/EntryGuidance.jl). 

The desired API is the following:

```julia
using CPEG

# initialize workspace
cp = cpeg_workspace()

# update atmosphere
update_atmosphere_params(cp,p::DensityParameters) 

# initialize trajectory 
cp.U = U # for some U 

# call call the alg
U_new = update_trajectory(cp,r,v) # cartesian position and velocity in MCMF
```
