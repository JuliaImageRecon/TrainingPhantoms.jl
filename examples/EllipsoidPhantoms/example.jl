using Pkg
Pkg.activate(joinpath(@__DIR__,".."))
Pkg.instantiate()
using TrainingPhantoms, Statistics
include(joinpath(@__DIR__,"..","helperFunctions","meanIntensityPlotAndIsosurface.jl"))
GLMakie.activate!()

# Generate a vessel phantom
rng = StableRNGs.StableRNG(1234)
N = (51,51,51)
phantom = ellipsoidPhantom(N; rng=rng)

# Visualize the phantom
meanIntensityPlotAndIsosurface(phantom)