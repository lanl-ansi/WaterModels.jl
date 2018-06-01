@testset "test minlp feasibility problems" begin
    #@testset "balerma" begin
    #    solution = run_feasibility("../test/data/epanet/balerma.inp", GenericWaterModel{StandardMINLPForm}, solver)
    #    @test solution["status"] == :LocalOptimal
    #end

    #@testset "foss_poly_1" begin
    #    solution = run_feasibility("../test/data/epanet/foss_poly_1.inp", GenericWaterModel{StandardMINLPForm}, solver)
    #    @test solution["status"] == :LocalOptimal
    #end

    #@testset "hanoi_extended" begin
    #    solution = run_feasibility("../test/data/epanet/hanoi_extended.inp", GenericWaterModel{StandardMINLPForm}, solver)
    #    @test solution["status"] == :LocalOptimal
    #end

    @testset "hanoi" begin
        solution = run_feasibility("../test/data/epanet/hanoi.inp", GenericWaterModel{StandardMINLPForm}, solver)
        @test solution["status"] == :LocalOptimal
    end

    #@testset "klmod" begin
    #    solution = run_feasibility("../test/data/epanet/klmod.inp", GenericWaterModel{StandardMINLPForm}, solver)
    #    @test solution["status"] == :LocalOptimal
    #end

    #@testset "rural" begin
    #    solution = run_feasibility("../test/data/epanet/rural.inp", GenericWaterModel{StandardMINLPForm}, solver)
    #    @test solution["status"] == :LocalOptimal
    #end

    #@testset "zj" begin
    #    solution = run_feasibility("../test/data/epanet/zj.inp", GenericWaterModel{StandardMINLPForm}, solver)
    #    @test solution["status"] == :LocalOptimal
    #end
end
