using Pkg, TrainingPhantoms, Statistics

# Install required packages
for P in [:CairoMakie, :GLMakie, :Makie, :StableRNGs]
    try
        @eval import $P
    catch
        Pkg.add(String(P))
        @eval import $P
    end
end
include(joinpath("..","helperFunctions","meanIntensityPlotAndIsosurface.jl"))
GLMakie.activate!()

# Generate a vessel phantom
rng = StableRNGs.StableRNG(1234)
N = (51,51,51)
phantom = ellipsoidPhantom(N; rng=rng)

# Visualize the phantom
meanIntensityPlotAndIsosurface(phantom)