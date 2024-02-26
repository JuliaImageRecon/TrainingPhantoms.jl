ellipsoidFunction(::NTuple{D, Int}) where D = throw(ArgumentError("Ellipsoid phantoms are currently not implemented for $D dimensions."))
ellipsoidFunction(::NTuple{2, Int}) = ellipse
ellipsoidFunction(::NTuple{3, Int}) = ellipsoid

# Based in ImagePhantoms.jl
function singleEllipsoid(N::NTuple{D,Int}, radius, shift, rotAngles) where D
  ranges = ntuple(d-> 1:N[d], D)
  objectFunction = ellipsoidFunction(N)
  ob = objectFunction(shift .+ (N.÷2), radius, rotAngles, 1.0f0)
  img = phantom(ranges..., [ob], 2)
  return img
end

"""
    scaleValue(x, a, b)
Scale a value `x` that is assumed to lie in the interval [0,1] 
to the interval [a,b] for given values a and b.
"""
function scaleValue(x::T, a::Real, b::Real=one(T)) where {T <: Real}
  @assert a >= zero(T) && a < one(T) "scaling factor a must lie in the interval [0,1) but is $a"
  @assert x >= zero(T) && x <= one(T) "value x is assumed to lie in the interval [0,1] but is $x"
  return a + (x*(b-a))
end

"""
  ellipsoidPhantom(
    N::NTuple{D,Int}; 
    rng::AbstractRNG = GLOBAL_RNG,
    numObjects::Int = rand(rng, 5:10),
    minRadius::Real=1.0,
    minValue::Real=0.1,
    allowOcclusion::Bool=false
  )

- `N`: size of the phantom image
- `rng`: random number generator
- `numObjects`: number of ellipses to generate
- `minRadius`: minimal radius in pixel
- `minValue`: minimal value of a single ellipse
- `allowOcclusion`: if `true` ellipse overshadows smaller values at its location, i.e., 
      new ellipses are not simply added to the exisiting object, instead the maximum is selected
"""
function ellipsoidPhantom(N::NTuple{D,Int}; rng::AbstractRNG = GLOBAL_RNG,
                          numObjects::Int = rand(rng, 5:10), minRadius::Real=1.0,
                          minValue::Real=0.1, allowOcclusion::Bool=false) where D
  img = zeros(N)
  @debug numObjects
  for m=1:numObjects
    # in 2D there is just one rotational degree of freedom
    rotAngles = ntuple(_ -> 2π*rand(rng), D == 2 ? 1 : D)
    radius = N .* scaleValue.(ntuple(_ -> 0.3*rand(rng), D), minRadius./N)
    shift = N .* ntuple(_ -> 0.6*(rand(rng)-0.5), D) 
    value = scaleValue.(rand(rng, 1), minValue)#.^2
    kernelWidth = ntuple(_ -> rand(rng)*N[1] / 20, D)
    P = singleEllipsoid(N, radius, shift, rotAngles)
    P = imfilter(P, Kernel.gaussian(kernelWidth))
    if allowOcclusion
      img = max.(img, value.*P)
    else
      img += value.*P
    end
  end
  return img
end