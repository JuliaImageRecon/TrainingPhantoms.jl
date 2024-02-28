# This file is originally based on Matlab code written by Christine Droigk

"""
vesselPath(N::NTuple{3,Int}; start, angle_xy, angle_xz, diameter, split_prob, change_prob, 
  max_change, splitnr, max_number_splits, stepsize, change_diameter_splitting, split_prob_factor, 
  change_prob_increase, rng)

Input parameters:
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
                             rng::AbstractRNG = GLOBAL_RNG)

  route = NTuple{3,Float64}[]
  push!(route, start)
  while all( 1 .<= route[end] .< N) # while the route of the vessel is inside the image area
    
    if (rand(rng, Float64) < change_prob) && (length(route) > 2)# if directional change of the vessel
        change_angle_xy = pi*(max_change-2*max_change*rand(rng,Float64)) # throw dice for the angles of directional change. For more or less variations replace (0.1-0.2*rand)
        change_angle_xz = pi*(max_change-2*max_change*rand(rng,Float64)) # same for other angle
        step_angle_xy = range(0, change_angle_xy, length=20) # to obtain a piecewise change of angle
        step_angle_xz = range(0, change_angle_xz, length=20)
        
        for i = 1:max(length(step_angle_xy),length(step_angle_xz)) # until the largest change of the two angles is reached
          if all( 1 .<= route[end] .< N) # check whether image boundaries are reached
            if i<=length(step_angle_xy) &&  i<=length(step_angle_xz) # both angles are still changing
              push!(route, route[end] .+ 
                    stepsize .* (cos(angle_xy + step_angle_xy[i])*cos(angle_xz + step_angle_xz[i]), 
                                 sin(angle_xy + step_angle_xy[i])*cos(angle_xz + step_angle_xz[i]), 
                                 sin(angle_xz + step_angle_xz[i])) )
            elseif i>length(step_angle_xy) &&  i<=length(step_angle_xz) # only xz-angle changes
              push!(route, route[end] .+ 
                    stepsize .* (cos(angle_xy + change_angle_xy)*cos(angle_xz + step_angle_xz[i]), 
                                 sin(angle_xy + change_angle_xy)*cos(angle_xz + step_angle_xz[i]), 
                                 sin(angle_xz + step_angle_xz[i])) )
            elseif i<=length(step_angle_xy) &&  i>length(step_angle_xz) # only xy-angle changes
              push!(route, route[end] .+ 
                    stepsize .* (cos(angle_xy+step_angle_xy[i])*cos(angle_xz + change_angle_xz), 
                                 sin(angle_xy+step_angle_xy[i])*cos(angle_xz + change_angle_xz), 
                                 sin(angle_xz + change_angle_xz)) )
            end
          end
        end
        # set current angles to new angles
        angle_xy = angle_xy + change_angle_xy # set current angles to new angles
        angle_xz = angle_xz + change_angle_xz
    else
        # if no directional change
        num_step = 15
        for _=1:num_step
          push!(route, route[end] .+ stepsize .* (cos(angle_xy)*cos(angle_xz), 
                                                sin(angle_xy)*cos(angle_xz), 
                                                sin(angle_xz))) # use old angle
        end
    end
    
    if rand(rng, Float64) < split_prob && (length(route) > 3) && (splitnr ≤ max_number_splits) # if vessel is splitting
        min_angle_diff = 0.15*pi
        angle_diff = rand(rng, Float64)*(pi/2 - min_angle_diff) + min_angle_diff # throw dice for angle between the resulting two vessel parts
        angle_diff_z = rand(rng, Float64)*(pi/2 - min_angle_diff) + min_angle_diff # same for z-angle 
      
        part_a = rand(rng, Float64) # Part of the angle related to the first division
        part_a_z = rand(rng, Float64) # Same for z-angle
        
        # Call recursively this function for each of the two vessel
        # segments with adjusted angles and diameter. Here, an adjustment
        # of the splitting probabilities and directional change
        # probabilities is made.
        routeA, diameterA = vesselPath(N; start=route[end], angle_xy=angle_xy-part_a*angle_diff, 
            angle_xz=angle_xz-part_a_z*angle_diff_z, diameter=diameter*change_diameter_splitting,
            split_prob=split_prob*split_prob_factor, change_prob=change_prob+change_prob_increase, max_change, splitnr=splitnr+1, 
            max_number_splits, stepsize, change_diameter_splitting, split_prob_factor, change_prob_increase, rng)
        routeB, diameterB = vesselPath(N; start=route[end], angle_xy=angle_xy+(1-part_a)*angle_diff,  
            angle_xz=angle_xy+(1-part_a_z)*angle_diff_z, diameter=diameter*change_diameter_splitting, 
            split_prob=split_prob*split_prob_factor, change_prob=change_prob+change_prob_increase, max_change, splitnr=splitnr+1, 
            max_number_splits, stepsize, change_diameter_splitting, split_prob_factor, change_prob_increase, rng)     
        
        if splitnr>1
          # set diameter for the splitting parts
          change_diameter_splitting_ = length(route) > 1 ? change_diameter_splitting : 1
          diameter_route = collect(range((1/change_diameter_splitting_)*diameter, diameter, length=length(route)))
        else
          diameter_route = diameter*ones(length(route))
        end
        append!(route, routeA, routeB)
        append!(diameter_route, diameterA, diameterB)
        return route, diameter_route
    end
  end

  if splitnr>1
    change_diameter_splitting_ = length(route) > 1 ? change_diameter_splitting : 1
    diameter_route = collect(range((1/change_diameter_splitting_)*diameter, diameter, length=length(route)))
  else
    diameter_route = diameter*ones(length(route))
  end

  return route, diameter_route
end


"""
vesselPhantom(N::NTuple{3,Int}; oversampling=2, kargs...)

Example usage:

  using GR, TrainingPhantoms, StableRNGs
  im = vesselPhantom((40,40,40); start=(1, 20, 20), angle_xy=0.0, angle_xz=0.0, 
                                diameter=2, split_prob=0.5, change_prob=0.5, max_change=0.2, splitnr=1,
                                rng = StableRNG(1));
  isosurface(im, isovalue=0.2, rotation=110, tilt=40)

"""
function vesselPhantom(N::NTuple{3,Int}; oversampling=2, rng = GLOBAL_RNG, kargs...)
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
