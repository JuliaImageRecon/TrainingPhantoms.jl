using Pkg
Pkg.activate(joinpath(@__DIR__,".."))
Pkg.instantiate()
using TrainingPhantoms, Statistics
include(joinpath(@__DIR__,"..","helperFunctions","meanIntensityPlotAndIsosurface.jl"))
GLMakie.activate!()

# Generate a vessel phantom
rng = StableRNGs.StableRNG(123)
N = (51,51,51)
phantom = vesselPhantom(N; 
    start=(1, N[2]/2, N[3]/2), angles=(0.0,0.0), 
    diameter=2.5, split_prob=0.4, change_prob=0.3, 
    max_change=0.3, splitnr=1, max_number_splits=1, rng=rng)

# Visualize the phantom
meanIntensityPlotAndIsosurface(phantom)