# This file is originally based on Matlab code written by Christine Droigk
"""
  appendRoute!(route, stepsize, angle_xy, angle_xz)

Append a new point to the route of the vessel.
"""
function appendRoute!(route, stepsize, angle_xy, angle_xz)
    push!(route, route[end] .+ stepsize .* (cos(angle_xy)*cos(angle_xz), 
                                            sin(angle_xy)*cos(angle_xz), 
                                            sin(angle_xz)))
end

"""
  changeDirection!(route, N, stepsize, angle_xy, angle_xz, change_prob, max_change, rng)

Simulate if and how the vessel changes its direction.
"""
function changeDirection!(route, N, stepsize, angle_xy, angle_xz, change_prob, max_change, steps_each_change, steps_no_change, rng)
  if (rand(rng, Float64) < change_prob) && (length(route) > 2)# if directional change of the vessel
    # randomly select new angles
    change_angle_xy = pi*(max_change-2*max_change*rand(rng,Float64))
    change_angle_xz = pi*(max_change-2*max_change*rand(rng,Float64))
    # create range such that change is gradually applied
    step_angle_xy = range(0, change_angle_xy, length=steps_each_change)
    step_angle_xz = range(0, change_angle_xz, length=steps_each_change)
    
    for i = 1:steps_each_change
      if all( 1 .<= route[end] .< N) # check whether image boundaries are reached
          appendRoute!(route, stepsize, angle_xy + step_angle_xy[i], angle_xz + step_angle_xz[i])
      end
    end
    # set current angles to new angles
    angle_xy = angle_xy + change_angle_xy
    angle_xz = angle_xz + change_angle_xz
  else
      # if no directional change
      for _=1:steps_no_change
        if all( 1 .<= route[end] .< N) # check whether image boundaries are reached
          appendRoute!(route, stepsize, angle_xy, angle_xz)          
        end
      end
  end
  return angle_xy, angle_xz
end

"""
  getDiameterRoute(route, diameter, change_diameter_splitting, splitnr)

Compute the diameter anlong the route of the vessel.
"""
function getDiameterRoute(route, diameter, change_diameter_splitting, splitnr)
  if splitnr>1
    # for the case that the vessel leaves the image during a split the change is set to 1
    # otherwise `range` throws an error since start and end do not match but length is 1
    change_diameter_splitting_ = length(route) > 1 ? change_diameter_splitting : 1
    # gradually change diameter from old value to current one
    diameter_route = collect(range((1/change_diameter_splitting_)*diameter, diameter, length=length(route)))
  else
    # no change if vessel did not already split once
    diameter_route = diameter*ones(length(route))
  end
  return diameter_route
end

"""
  vesselPath(N::NTuple{3,Int}; start, angle_xy, angle_xz, diameter, split_prob, change_prob, 
    max_change, splitnr, max_number_splits, stepsize, change_diameter_splitting, split_prob_factor, 
    change_prob_increase, rng)

### Input parameters:
* N: Image size, given as a 3 tuple
* start: starting point given as a 3x1 vector
* angle_xy: angle in radians describing the starting angle with which the
* vessel runs. 
* angle_xz: same, but angle from xy-plane to vessel's z
* diameter: starting diameter of vessel
* split_prob: probability for a splitting of the vessel into two vessel segments. Values between 0 and 1.
* change_prob: probability for directional change of the vessel route. Values between 0 and 1.
* max_change: max_change * pi specifies the maximum direction-change angle.
* splitnr: used for recursive call of the function. For the first call set it to 1. 
* max_number_splits: maximum number of splits of the vessel.
* stepsize: stepsize of the vessel.
* change_diameter_splitting: Indicates by how much the diameter decreases when the vessel splits
* split_prob_factor: Factor by which the split probability `split_prob` is multiplied when the vessel splits
* change_prob_increase: Increase of the change probability `change_prob` when the vessel splits
* steps_each_change: Number of steps for the change of the vessel direction
* steps_no_change: Number of steps for the case that the vessel does not change its direction
* rng: Random number generator
 
Output:
* route: A length N vector containing the points 3D of the route of the vessel. The
* length N depends on the random route.
* diameter_route: A length N vector containing the diameter of the vessel at the positions of the route.
"""
function vesselPath(N::NTuple{3,Int}; 
                             start, 
                             angle_xy, 
                             angle_xz, 
                             diameter, 
                             split_prob, 
                             change_prob, 
                             max_change, 
                             splitnr,
                             max_number_splits = Inf,
                             stepsize = max(N...)/200,
                             change_diameter_splitting = 0.6,
                             split_prob_factor = 0.5,
                             change_prob_increase = 0.01,
                             steps_each_change = 20,
                             steps_no_change = 15,
                             rng::AbstractRNG = GLOBAL_RNG)

  route = NTuple{3,Float64}[]
  push!(route, start)
  while all( 1 .<= route[end] .< N) # while the route of the vessel is inside the image area
    
    angle_xy, angle_xz =  changeDirection!(route, N, stepsize, angle_xy, angle_xz, change_prob, 
        max_change, steps_each_change, steps_no_change, rng)
    
    # Splitting part
    if rand(rng, Float64) < split_prob && (length(route) > 3) && (splitnr ≤ max_number_splits) # if vessel is splitting
        min_angle_diff = 0.15*pi
        # Randomly select the angle between the resulting two vessel parts
        angle_diff_xy = rand(rng, Float64)*(pi/2 - min_angle_diff) + min_angle_diff
        angle_diff_xz = rand(rng, Float64)*(pi/2 - min_angle_diff) + min_angle_diff
        # Randomly select the part of the angle related to the first division
        part_a_xy = rand(rng, Float64)
        part_a_xz = rand(rng, Float64)
        
        # Call recursively this function for each of the two vessel
        # segments with adjusted angles and diameter. Here, an adjustment
        # of the splitting probabilities and directional change
        # probabilities is made.
        routeA, diameterA = vesselPath(N; start=route[end], angle_xy=angle_xy-part_a_xy*angle_diff_xy, 
            angle_xz=angle_xz-part_a_xz*angle_diff_xz, diameter=diameter*change_diameter_splitting,
            split_prob=split_prob*split_prob_factor, change_prob=change_prob+change_prob_increase, 
            max_change, splitnr=splitnr+1, max_number_splits, stepsize, change_diameter_splitting, 
            split_prob_factor, change_prob_increase, steps_each_change, steps_no_change, rng)
        routeB, diameterB = vesselPath(N; start=route[end], angle_xy=angle_xy+(1-part_a_xy)*angle_diff_xy,  
            angle_xz=angle_xy+(1-part_a_xz)*angle_diff_xz, diameter=diameter*change_diameter_splitting, 
            split_prob=split_prob*split_prob_factor, change_prob=change_prob+change_prob_increase, 
            max_change, splitnr=splitnr+1, max_number_splits, stepsize, change_diameter_splitting, 
            split_prob_factor, change_prob_increase, steps_each_change, steps_no_change, rng)     
        
        diameter_route = getDiameterRoute(route, diameter, change_diameter_splitting, splitnr)
        append!(route, routeA, routeB)
        append!(diameter_route, diameterA, diameterB)
        return route, diameter_route
    end
  end
  diameter_route = getDiameterRoute(route, diameter, change_diameter_splitting, splitnr)
  return route, diameter_route
end


"""
  vesselPhantom(N::NTuple{3,Int}; oversampling=2, kargs...)

### Input parameters:
* N: Image size, given as a 3 tuple
* oversampling: Oversampling factor for the phantom. Default is 2.
* rng: Random number generator
* kernelWidth: Width (standard deviation) of the Gaussian kernel (in pixel) used for smoothing the phantom. If nothing is given, a random value is chosen.
* kargs...: remaining keyword arguments for `vesselPath`

### Example usage:

  using GLMakie, TrainingPhantoms, StableRNGs

  im = vesselPhantom((51,51,51); 
          start=(1, 25, 25), angle_xy=0.0, angle_xz=0.0, 
          diameter=2.5, split_prob=0.4, change_prob=0.3, 
          max_change=0.3, splitnr=1, max_number_splits=1, rng StableRNG(123));

  f = Figure(size=(300,300))
  ax = Axis3(f[1,1], aspect=:data)
  volume!(ax, im, algorithm=:iso, isorange=0.13, isovalue=0.3, colormap=:viridis, colorrange=[0.0,0.2])
"""
function vesselPhantom(N::NTuple{3,Int}; oversampling=2, rng = GLOBAL_RNG, kernelWidth=nothing, kargs...)
  route, diameter_route = vesselPath(N; rng, kargs...)
  # add small sphere for every entry in the route
  obs = [ sphere( Float32.(route[i]), Float32(diameter_route[i]), 1.0f0) for i=eachindex(route) ]
  ranges = ntuple(d-> 1:N[d], 3)
  img = phantom(ranges..., obs, oversampling)
  if isnothing(kernelWidth)
    # filter for smoothing, offset to ensure minimal filter width
    filterWidth = (1.0-0.3)*rand(rng) + 0.3
    kernelWidth = ntuple(_ -> filterWidth*N[1] / 20, 3)
  end  
  img = imfilter(img, Kernel.gaussian(kernelWidth))
  img[img .> 1] .= 1
  return img
end
