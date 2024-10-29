# This file is originally based on Matlab code written by Christine Droigk
"""
    appendRoute!(route, stepsize, angles)

Append a new point to the route of the vessel.
"""
function appendRoute!(route, stepsize, angles::NTuple{2,Float64})
  angle_xy, angle_xz = angles
  push!(route, route[end] .+ stepsize .* (cos(angle_xy)*cos(angle_xz), 
                                          sin(angle_xy)*cos(angle_xz), 
                                          sin(angle_xz)))
end

function appendRoute!(route, stepsize, angles::NTuple{1,Float64})
  (angle_xy,) = angles
  push!(route, route[end] .+ stepsize .* (cos(angle_xy), 
                                          sin(angle_xy)))
end

"""
    changeDirection!(route, N, stepsize, angles, changeProb, maxChange, rng)

Simulate if and how the vessel changes its direction.
"""
function changeDirection!(N::NTuple{D,Int}, route, stepsize, angles, changeProb, maxChange, stepsEachChange, stepsNoChange, rng) where D
  if (rand(rng, Float64) < changeProb) && (length(route) > 2)# if directional change of the vessel
    # randomly select new angles
    change_angles = ntuple(_ -> 2*maxChange*rand(rng,Float64), D-1)
    # create range such that change is gradually applied
    step_angles = zip(ntuple(d -> range(0, change_angles[d], length=stepsEachChange), D-1)...) |> collect
    for i = 1:stepsEachChange
      if all( 1 .<= route[end] .<= N) # check whether image boundaries are reached
          appendRoute!(route, stepsize, angles .+ step_angles[i])
      end
    end
    # set current angles to new angles
    angles = angles .+ change_angles
  else
      # if no directional change
      for _=1:stepsNoChange
        if all( 1 .<= route[end] .<= N) # check whether image boundaries are reached
          appendRoute!(route, stepsize, angles)          
        end
      end
  end
  return angles
end

"""
    getDiameterRoute(route, diameter, splitDiameterChange, splitnr)

Compute the diameter anlong the route of the vessel.
"""
function getDiameterRoute(route, diameter, splitDiameterChange, splitnr)
  if splitnr>1
    # for the case that the vessel leaves the image during a split the change is set to 1
    # otherwise `range` throws an error since start and end do not match but length is 1
    splitDiameterChange_ = length(route) > 1 ? splitDiameterChange : 1
    # gradually change diameter from old value to current one
    diameterRoute = collect(range((1/splitDiameterChange_)*diameter, diameter, length=length(route)))
  else
    # no change if vessel did not already split once
    diameterRoute = diameter*ones(length(route))
  end
  return diameterRoute
end

"""
    vesselPath(N::NTuple{D,Int}; start, angles, diameter, splitProb, changeProb, 
      maxChange, splitnr, maxNumSplits, stepsize, splitDiameterChange, splitProbFactor, 
      changeProbIncrease, stepsEachChange, stepsNoChange, rng)


# Arguments
- `N`: Image size, given as a D tuple
- `start`: starting point given as a D tuple
- `angles`: D-1 tuple of angles for the vessel's direction (xy in 2D, xy and xz in 3D)
- `diameter`: starting diameter of vessel
- `splitProb`: probability for a splitting of the vessel into two vessel segments. Values between 0 and 1.
- `changeProb`: probability for directional change of the vessel route. Values between 0 and 1.
- `maxChange`: maxChange * pi specifies the maximum direction-change angle.
- `splitnr`: used for recursive call of the function. For the first call set it to 1. 
- `maxNumSplits`: maximum number of splits of the vessel.
- `stepsize`: stepsize of the vessel.
- `splitDiameterChange`: Indicates by how much the diameter decreases when the vessel splits
- `splitProbFactor`: Factor by which the split probability `splitProb` is multiplied when the vessel splits
- `changeProbIncrease`: Increase of the change probability `changeProb` when the vessel splits
- `stepsEachChange`: Number of steps for the change of the vessel direction
- `stepsNoChange`: Number of steps for the case that the vessel does not change its direction
- `rng`: Random number generator
 
# Output
- `route`: A length N vector containing the D dimensional points of the route of the vessel. The 
length N depends on the random route.
- `diameterRoute`: A length N vector containing the diameter of the vessel at the positions of the route.
"""
function vesselPath(N::NTuple{D,Int}; 
                    start,
                    angles::NTuple{Da,Real},
                    diameter::Real, 
                    splitProb::Real, 
                    changeProb::Real, 
                    maxChange::Real, 
                    splitnr::Int=1,
                    maxNumSplits::Real = Inf,
                    stepsize::Real = max(N...)/200,
                    splitDiameterChange::Real = 0.6,
                    splitProbFactor::Real = 0.5,
                    changeProbIncrease::Real = 0.01,
                    stepsEachChange::Int = 20,
                    stepsNoChange::Int = 15,
                    rng::AbstractRNG = GLOBAL_RNG) where {D, Da}
  @assert (D-1) == Da "The length of the angles tuple must be $(D-1) but is $Da"
  route = NTuple{D,Float64}[]
  push!(route, start)
  while all( 1 .<= route[end] .<= N) # while the route of the vessel is inside the image area
    angles =  changeDirection!(N, route, stepsize, angles, changeProb, 
        maxChange, stepsEachChange, stepsNoChange, rng)
    
    # Splitting part
    if rand(rng, Float64) < splitProb && (length(route) > 3) && (splitnr ≤ maxNumSplits) # if vessel is splitting
        min_angle_diff = 0.15*pi
        # Randomly select the angle between the resulting two vessel parts
        angle_diffs = ntuple(_ -> rand(rng, Float64)*(pi/2 - min_angle_diff) + min_angle_diff, D-1)
        # Randomly select the part of the angle related to the first division
        parts_a = ntuple(_ -> rand(rng, Float64), D-1)
        
        # Call recursively this function for each of the two vessel
        # segments with adjusted angles and diameter. Here, an adjustment
        # of the splitting probabilities and directional change
        # probabilities is made.    
        routeA, diameterA = vesselPath(N; start=route[end], angles=angles .- parts_a .* angle_diffs,
            diameter=diameter*splitDiameterChange,
            splitProb=splitProb*splitProbFactor, changeProb=changeProb+changeProbIncrease, 
            maxChange, splitnr=splitnr+1, maxNumSplits, stepsize, splitDiameterChange, 
            splitProbFactor, changeProbIncrease, stepsEachChange, stepsNoChange, rng)
        routeB, diameterB = vesselPath(N; start=route[end], angles=angles .+ (1 .- parts_a) .* angle_diffs,  
            diameter=diameter*splitDiameterChange, 
            splitProb=splitProb*splitProbFactor, changeProb=changeProb+changeProbIncrease, 
            maxChange, splitnr=splitnr+1, maxNumSplits, stepsize, splitDiameterChange, 
            splitProbFactor, changeProbIncrease, stepsEachChange, stepsNoChange, rng)     
        
        diameterRoute = getDiameterRoute(route, diameter, splitDiameterChange, splitnr)
        append!(route, routeA, routeB)
        append!(diameterRoute, diameterA, diameterB)
        return route, diameterRoute
    end
  end
  diameterRoute = getDiameterRoute(route, diameter, splitDiameterChange, splitnr)
  return route, diameterRoute
end

sphereFunction(::NTuple{D, Int}) where D = throw(ArgumentError("Vessel phantoms are currently not implemented for $D dimensions."))
sphereFunction(::NTuple{2, Int}) = circle
sphereFunction(::NTuple{3, Int}) = sphere

"""
    vesselPhantom(N::NTuple{D,Int}; oversampling=2, kargs...)

Generate a random vessel phantom.

# Arguments
* N: Image size, given as a D tuple
* oversampling: Oversampling factor for the phantom. Default is 2.
* rng: Random number generator
* kernelWidth: Width (standard deviation) of the Gaussian kernel (in pixel) used for smoothing the phantom. If nothing is given, a random value is chosen.
* kargs...: remaining keyword arguments for `vesselPath`

# Examples

  using GLMakie, TrainingPhantoms, StableRNGs

  im = vesselPhantom((51,51,51); 
          start=(1, 25, 25), angles=(0.0, 0.0), 
          diameter=2.5, splitProb=0.4, changeProb=0.3, 
          maxChange=0.3, splitnr=1, maxNumSplits=1, rng StableRNG(123));

  f = Figure(size=(300,300))
  ax = Axis3(f[1,1], aspect=:data)
  volume!(ax, im, algorithm=:iso, isorange=0.13, isovalue=0.3, colormap=:viridis, colorrange=[0.0,0.2])
"""
function vesselPhantom(N::NTuple{D,Int}; oversampling=2, rng = GLOBAL_RNG, kernelWidth=nothing, kargs...) where D
  objectFunction = sphereFunction(N)
  route, diameterRoute = vesselPath(N; rng, kargs...)
  # add small sphere for every entry in the route
  obs = [ objectFunction( Float32.(route[i]), Float32(diameterRoute[i]), 1.0f0) for i=eachindex(route) ]
  ranges = ntuple(d-> 1:N[d], D)
  img = phantom(ranges..., obs, oversampling)
  if isnothing(kernelWidth)
    # filter for smoothing, offset to ensure minimal filter width
    filterWidth = (1.0-0.3)*rand(rng) + 0.3
    kernelWidth = ntuple(_ -> filterWidth*N[1] / 20, D)
  end  
  img = imfilter(img, Kernel.gaussian(kernelWidth))
  img[img .> 1] .= 1
  img[img .< 0] .= 0
  return img
end
