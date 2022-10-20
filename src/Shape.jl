
#= Based in ImagePhantoms.jl
function makeEllipse(N::NTuple{D,Int}, radius, shift, rot) where D

  ranges = ntuple(d-> 1:N[d], D)
  ob = ellipsoid(Float32.(Tuple(shift .+ (N.÷2)) ), Float32.(Tuple(radius)), (Float32(rand(Float32)*2*pi),Float32(rand(Float32)*2*pi*0) ), 1.0f0)
  img = phantom(ranges..., [ob], 2)
  return img
end
=#

function singleEllipsoid(N::NTuple{D,Int}, radius, shift, rot) where D

  I = zeros(N)
  N_ = collect(Tuple(N))
  R = N_ .* radius

  for i in CartesianIndices(N)
    r = collect(Tuple(i)) 
    I[i] = norm(  (rot*(r .- (N.÷2)) ./ R) .- shift ./radius)  <= 1.0
  end

  return I
end

function ellipsoidPhantom(N::NTuple{D,Int}; rng::AbstractRNG = GLOBAL_RNG,
                                            numObjects::Int = rand(rng, 5:10)) where D
  img = zeros(N)
  for m=1:numObjects
    rot = rand(rng, RotMatrix{D})
    radius = rand(rng, D)*0.3
    shift = rand(rng, D)*0.3 
    value = (rand(rng, 1)).^2
    kernelWidth = ntuple(d-> rand(rng)*N[1] / 20, D)
    if shift > radius
      shift = radius
    end
    P = singleEllipsoid(N, radius, shift, rot)
    P = imfilter(P, Kernel.gaussian(kernelWidth))
    img += value.*P
  end
  return img
end