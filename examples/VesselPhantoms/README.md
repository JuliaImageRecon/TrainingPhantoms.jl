# Vessel Phantom Example
This example generates a vessel phantom and visualizes it using a mean intensity plot and an isosurface volume rendering.

## Installation
To use this code, you must first install `TrainingPhantoms` by running
```julia
using Pkg
Pkg.add("TrainingPhantoms")
```
Afterwards, you can load the package by entering
```julia
using TrainingPhantoms
```

## Execution
Finally, just include the example code
```julia
include(joinpath(dirname(pathof(TrainingPhantoms)),"..","VesselPhantoms", "vesselPhantom.jl"))
```
After the vessel was generated, the phantom will be shown in a new window.