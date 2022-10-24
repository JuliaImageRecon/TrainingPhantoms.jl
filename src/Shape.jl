
# Based in ImagePhantoms.jl
function singleEllipsoid(N::NTuple{D,Int}, radius, shift, rotAngles) where D

  ranges = ntuple(d-> 1:N[d], D)
  ob = ellipsoid(shift .+ (N.÷2), radius, rotAngles, 1.0f0)
  img = phantom(ranges..., [ob], 2)
  return img
end

function ellipsoidPhantom(N::NTuple{D,Int}; rng::AbstractRNG = GLOBAL_RNG,
                                            numObjects::Int = rand(rng, 5:10)) where D
  img = zeros(N)
  @info numObjects
  for m=1:numObjects
    rotAngles = ntuple(_ -> 2π*rand(rng), D)
    radius = N .* ntuple(_ -> 0.3*rand(rng), D)
    shift = N .* ntuple(_ -> 0.6*(rand(rng)-0.5), D) 
    value = (rand(rng, 1))#.^2
    kernelWidth = ntuple(_ -> rand(rng)*N[1] / 20, D)
    if shift > radius
      shift = radius
    end
    P = singleEllipsoid(N, radius, shift, rotAngles)
    P = imfilter(P, Kernel.gaussian(kernelWidth))
    img += value.*P
  end
  return img
end