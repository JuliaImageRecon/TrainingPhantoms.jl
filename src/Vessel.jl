# This file is originally based on Matlab code written by Christine Droigk

"""
vesselPath(N::NTuple{3,Int}; start, angle_xy, angle_xz, diameter, split_prob, change_prob, max_change, splitnr)

Input parameters:
* start: starting point given as a 3x1 vector
* angle_xy: angle in radians describing the starting angle with which the
* vessel runs. 
* angle_xz: same, but angle from xy-plane to vessel's z
* diameter: starting diameter of vessel
* split_prob: probability for a splitting of the vessel into two vessel segments. Values between 0 and 1.
* change_prob: probability for directional change of the vessel route. Values between 0 and 1.
* max_change: max_change * pi specifies the maximum direction-change angle.
* splitnr: used for recursive call of the function. For the first call set it to 1. 
* N: Image size, given as a 3x1 vector
 
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
                             rng::AbstractRNG = GLOBAL_RNG)

  stepsize = 0.25
  change_diameter_splitting = 4/5 # Indicates by how much the diameter decreases when the vessel is divided
  route = NTuple{3,Float64}[]
  push!(route, start)

  counter = 1;
  while all( 1 .<= route[end] .< N) # while the route of the vessel is inside the image area
    
    if rand(rng, Float64) < change_prob # if directional change of the vessel
        change_angle_xy = pi*(max_change-2*max_change*rand(rng,Float64)) # throw dice for the angles of directional change. For more or less variations replace (0.1-0.2*rand)
        change_angle_xz = pi*(max_change-2*max_change*rand(rng,Float64)) # same for other angle
        step_angle_xy = 0:(sign(change_angle_xy)*0.1*stepsize):change_angle_xy # to obtain a piecewise change of angle
        step_angle_xz = 0:(sign(change_angle_xz)*0.1*stepsize):change_angle_xz
        
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
        push!(route, route[end] .+ stepsize .* (cos(angle_xy)*cos(angle_xz), 
                                                sin(angle_xy)*cos(angle_xz), 
                                                sin(angle_xz))) # use old angle
    end
    
    
    if rand(rng, Float64) < split_prob && counter > 5 # if vessel is splitting
        angle_diff = rand(rng, Float64)*pi/2 # throw dice for angle between the resulting two vessel parts
        angle_diff_z = rand(rng, Float64)*pi/2 # same for z-angle 
      
        part_a = rand(rng, Float64) # Part of the angle related to the first division
        part_a_z = rand(rng, Float64) # Same for z-angle
        
        # Call recursively this function for each of the two vessel
        # segments with adjusted angles and diameter. Here, an adjustment
        # of the splitting probabilities and directional change
        # probabilities is made.
        routeA, diameterA = vesselPath(N; start=route[end], angle_xy=angle_xy-part_a*angle_diff, 
            angle_xz=angle_xz-part_a_z*angle_diff_z, diameter=diameter*change_diameter_splitting,
            split_prob=split_prob/2, change_prob=change_prob+0.01, max_change, splitnr=splitnr+1, rng)
        routeB, diameterB = vesselPath(N; start=route[end], angle_xy=angle_xy+(1-part_a)*angle_diff,  
            angle_xz=angle_xy+(1-part_a_z)*angle_diff_z, diameter=diameter*change_diameter_splitting, 
            split_prob=split_prob/2, change_prob=change_prob+0.01, max_change, splitnr=splitnr+1, rng)     
        
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
    
    counter = counter+1
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
function vesselPhantom(N::NTuple{3,Int}; oversampling=2, kargs...)
  Q = zeros(N);

  route, diameter_route = vesselPath(N; kargs...)

  obs = [ ellipsoid( Float32.(route[i]), Float32.(ntuple(_->diameter_route[i],3)), 
                              (0,0), 1.0f0) for i=1:length(route)]
  ranges = ntuple(d-> 1:N[d], 3)
  img = phantom(ranges..., obs, oversampling)
  img[img .> 1] .= 1
  return img
end
