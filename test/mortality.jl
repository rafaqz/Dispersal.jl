using Dispersal
using Test
using DynamicGrids

@testset "Mortality without resistance grid" begin
    init = (
        exposure   = [100.0 0.0;
                      100.0 0.0],
        population = [100.0 100.0;
                      0.0   0.0],
    )

    @testset "LoglogisticMortality" begin
        loglogmortalityrule = LoglogisticMortality{Tuple{:population,:exposure},:population}(
            median=20, hillcoefficient=0.5, timestep=1
        )
        output = ArrayOutput(init; tspan=1:3)
        sim!(output, loglogmortalityrule)
        @test output[1].exposure == output[2].exposure == output[3].exposure
        @test output[2].population ≈ [30.9017 100.0; 0.0 0.0] atol=1e-4
        @test output[3].population ≈ [9.54915 100.0; 0.0 0.0] atol=1e-4
    end

    @testset "ExponentialMortality low threshold" begin
        expmortalityrule = ExponentialMortality{Tuple{:population,:exposure},:population}(
            rate=0.01, threshold = 20.0, timestep=1
        )
        output = ArrayOutput(init; tspan=1:3)
        sim!(output, expmortalityrule)
        @test output[1].exposure == output[2].exposure == output[3].exposure
        @test output[2].population ≈ [44.9329 100.0; 0.0 0.0] atol=1e-4
        @test output[3].population ≈ [20.1897 100.0; 0.0 0.0] atol=1e-4
    end
    @testset "ExponentialMortality high threshold" begin
        expmortalityrule = ExponentialMortality{Tuple{:population,:exposure},:population}(
            rate=0.01, threshold=200.0, timestep=1
        )
        output = ArrayOutput(init; tspan=1:3)
        sim!(output, expmortalityrule)
        @test output[1].exposure == output[2].exposure == output[3].exposure
        @test output[1].population ==  output[2].population ==  output[3].population
    end
end

@testset "Mortality with Resistance grid" begin
    
    phenotyperule = Cell{:phenotype, :phenotype}() do data, state, I
        return state * 10
    end 
    
    @testset "LoglogisticMortality with Grid median" begin

        init = (
            exposure   = [100.0 0.0;
                          100.0 0.0],
            population = [100.0 100.0;
                          0.0   0.0],
            phenotype  = [20 20;
                          20 20], 
        )
        loglogmortalityrule = LoglogisticMortality{Tuple{:population,:exposure},:population}(
            median=Grid(:phenotype), hillcoefficient=0.5, timestep=1
        )
        output = ArrayOutput(init; tspan=1:3)
        sim!(output, Ruleset(loglogmortalityrule, phenotyperule))
        @test output[1].exposure == output[2].exposure == output[3].exposure
        @test output[2].population ≈ [30.9016 100.0; 0.0 0.0] atol=1e-4
        @test output[3].population ≈ [18.1017 100.0; 0.0 0.0] atol=1e-4
    end

    @testset "ExponentialMortality with Grid rate" begin
        init = (
            exposure   = [100.0 0.0;
                          100.0 0.0],
            population = [100.0 100.0;
                          0.0   0.0],
            phenotype  = [0.01 0.01;
                          0.01 0.01], 
        )
 
        expmortalityrule = ExponentialMortality{Tuple{:population,:exposure},:population}(
            rate=Grid(:phenotype), threshold=20.0, timestep=1
        )
        output = ArrayOutput(init; tspan=1:3)
        sim!(output, Ruleset(expmortalityrule, phenotyperule))
        @test output[1].exposure == output[2].exposure == output[3].exposure
        @test output[2].population ≈ [44.9328 100.0; 0.0 0.0] atol=1e-4
        @test output[3].population ≈ [0.01507 100.0; 0.0 0.0] atol=1e-4
    end
end

