using TrainingPhantoms
using StableRNGs
using Test

@testset "TrainingPhantoms.jl" begin
  N = (40,40,40)
  @testset "Vessel" begin    
    @testset "getDiameterRoute" begin
      route = [1, 2, 3, 4, 5]
      diameter = 2.0
      splitDiameterChange = 0.5
      splitnr = 2
  
      @testset "Normal operation" begin
          result = TrainingPhantoms.getDiameterRoute(route, diameter, splitDiameterChange, splitnr)
          @test length(result) == length(route)
          @test result[1] == (1/splitDiameterChange)*diameter
          @test result[end] == diameter
      end
  
      @testset "No split" begin
          result = TrainingPhantoms.getDiameterRoute(route, diameter, splitDiameterChange, 1)
          @test all(result .== diameter)
      end
    end
  end

  @testset "vesselPath Tests" begin
    rng = StableRNG(1)
    @testset "Normal operation" begin
        route, diameterRoute = TrainingPhantoms.vesselPath(N; start=(1,20,20), angles=(0.0,0.0), diameter=2, 
                                           splitProb=0.5, changeProb=0.5, maxChange=0.2, 
                                           splitnr=1, maxNumSplits=2, stepsize=0.25, 
                                           splitDiameterChange=4/5, splitProbFactor=0.5, 
                                           changeProbIncrease=0.01, rng=rng)
        @test length(route) == length(diameterRoute)
    end
    @testset "Start outside volume" begin
        route, diameterRoute = TrainingPhantoms.vesselPath(N; start=(1,50,20), angles=(0.0,0.0), diameter=2, 
                                           splitProb=0.5, changeProb=0.5, maxChange=0.2, 
                                           splitnr=1, maxNumSplits=2, stepsize=0.25, 
                                           splitDiameterChange=4/5, splitProbFactor=0.5, 
                                           changeProbIncrease=0.01, rng=rng)
        @test length(route) == 1
        @test length(diameterRoute) == 1
    end
  end
  @testset "vesselPhantom" begin
    im = vesselPhantom(N; start=(1, 20, 20), angles = (0.0, 0.0),
                       diameter=2, splitProb=0.5, changeProb=0.5, maxChange=0.2, splitnr=1,
                       rng = StableRNG(1));

    im2 = vesselPhantom(N; start=(1, 20, 20), angles = (0.0, 0.0), 
                        diameter=2, splitProb=0.5, changeProb=0.5, maxChange=0.2, splitnr=1,
                        rng = StableRNG(1));
    @test im ≈ im2
    @test size(im) == N
    @test maximum(im) == 1.0
    @test minimum(im) == 0.0

    im = vesselPhantom((40,40); start=(1, 20), angles = (0.0,), 
                       diameter=2, splitProb=0.5, changeProb=0.5, maxChange=0.2, splitnr=1,
                       rng = StableRNG(1));
    @test size(im) == (40,40)

    @test_throws ArgumentError vesselPhantom((20,20,20,20))

  end

  @testset "Ellipsoid" begin
    # Ellipsoid Phantom
    im = ellipsoidPhantom(N; rng=StableRNG(1))
    im2 = ellipsoidPhantom(N; rng=StableRNG(1))
    @test im ≈ im2
    @test size(im) == N

    im = ellipsoidPhantom((20,20); allowOcclusion=true)
    @test maximum(im) <= 1

    @test_throws ArgumentError ellipsoidPhantom((20,20,20,20))
    @test_throws ArgumentError TrainingPhantoms.ellipsoidFunction((20,20,20,20))
  end
end
