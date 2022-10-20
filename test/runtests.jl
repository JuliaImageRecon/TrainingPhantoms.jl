using TrainingPhantoms
using StableRNGs
using Test

@testset "TrainingPhantoms.jl" begin
  N = (40,40,40)

  # Vessel Phantom
  im = vesselPhantom(N; start=(1, 20, 20), angle_xy=0.0, angle_xz=0.0, 
                                diameter=2, split_prob=0.5, change_prob=0.5, max_change=0.2, splitnr=1,
                                rng = StableRNG(1));

  im2 = vesselPhantom(N; start=(1, 20, 20), angle_xy=0.0, angle_xz=0.0, 
                                diameter=2, split_prob=0.5, change_prob=0.5, max_change=0.2, splitnr=1,
                                rng = StableRNG(1));
  @test im ≈ im2

   # Ellipsoid Phantom
  im = ellipsoidPhantom(N; rng=StableRNG(1))
  im2 = ellipsoidPhantom(N; rng=StableRNG(1))
  @test im ≈ im2

  #isosurface(im, isovalue=0.2, rotation=110)
end
