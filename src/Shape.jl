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
    getRotationMatrix(D::Int, rotAngles::NTuple{Dr, <:Real}) where Dr
Get the rotation matrix for the given rotation angles.

### Input parameters:
- `D`: dimension of the object
- `rotAngles`: rotation angles
"""
function getRotationMatrix(D::Int, rotAngles::NTuple{Dr, <:Real}) where Dr
  R = zeros(Float64, (D,D))
  for i=1:D
    # readout columns of rotation matrix with unit vectors
    R[:,i] = [rotate_vector(ntuple(j -> (i==j ? 1.0 : 0.0), D), rotAngles)...]
  end
  return R
end

"""
    getEllipsoidMatrix(radius::NTuple{Dr, <:Real}, rotAngles::NTuple{Da, <:Real}) where {Dr, Da}
Get the matrix that describes the ellipsoid in the rotated coordinate system.

### Input parameters:
- `radius`: radii of the ellipsoid
- `rotAngles`: rotation angles
"""
function getEllipsoidMatrix(radius::NTuple{Dr, <:Real}, rotAngles::NTuple{Da, <:Real}) where {Dr, Da}
  R = getRotationMatrix(Dr, rotAngles)
  A = Diagonal([1/radius[i]^2 for i=1:Dr])
  return transpose(R)*A*R
end

"""
    ellipsoidBoundingBox(radius::NTuple{Dr, <:Real}, rotAngles::NTuple{Da, <:Real}) where {Dr, Da}
Get the bounding box of the specified ellipsoid.

### Input parameters:
- `radius`: radii of the ellipsoid
- `rotAngles`: rotation angles
"""
function ellipsoidBoundingBox(radius::NTuple{Dr, <:Real}, rotAngles::NTuple{Da, <:Real}) where {Dr, Da}
  B = getEllipsoidMatrix(radius, rotAngles) |> inv
  return sqrt.(diag(B))
end

# rotate vector v by the given rotation angles
rotate_vector(v::NTuple{3, <:Real}, rotAngles::NTuple{3, <:Real}) = ImagePhantoms.Rxyz_inv(v, rotAngles...)
rotate_vector(v::NTuple{2, <:Real}, rotAngles::NTuple{1, <:Real}) = ImagePhantoms.rotate2d(v, rotAngles...)

# Get the number of rotational degrees of freedom (according to Euclidean group)
rotDOF(::NTuple{D,<:Real}) where D = Int(D*(D-1)/2)

"""
  ellipsoidPhantom(
    N::NTuple{D,Int}; 
    rng::AbstractRNG = GLOBAL_RNG,
    numObjects::Int = rand(rng, 5:10),
    minRadius::Real=1.0,
    minValue::Real=0.1,
    allowOcclusion::Bool=false
  )

### Input parameters:
- `N`: size of the phantom image
- `rng`: random number generator
- `numObjects`: number of ellipses to generate
- `minRadius`: minimal radius in pixel
- `minValue`: minimal value of a single ellipse
- `allowOcclusion`: if `true` ellipse overshadows smaller values at its location, i.e., 
      new ellipses are not simply added to the exisiting object, instead the maximum is selected
- `pixelMargin`: minimal distance of the object to the edge of the image
"""
function ellipsoidPhantom(N::NTuple{D,Int}; rng::AbstractRNG = GLOBAL_RNG,
                          numObjects::Int = rand(rng, 5:10), minRadiusPixel::Real=1.0,
                          maxRadiusPercent::Real=0.3, maxShiftPercent::Real=0.6,
                          minValue::Real=0.1, allowOcclusion::Bool=false, pixelMargin::Real=1) where D
  img = zeros(N)
  @debug numObjects
  for m=1:numObjects
    rotAngles = ntuple(_ -> 2π*rand(rng), rotDOF(N))
    radius = N .* scaleValue.(ntuple(_ -> rand(rng), D), minRadiusPixel./N, maxRadiusPercent)
    shift = N .* ntuple(_ -> maxShiftPercent*(rand(rng)-0.5), D)
    boundingBox = ellipsoidBoundingBox(radius, rotAngles)
    # if the bounding box is too close to the edge of the image, shift the object to have `pixelMargin` distance
    shift = ntuple(i -> (boundingBox[i]+abs(shift[i]) < (N[i]/2 - pixelMargin)) ? shift[i] : (N[i]/2 - boundingBox[i] - pixelMargin)*sign(shift[i]), D)
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