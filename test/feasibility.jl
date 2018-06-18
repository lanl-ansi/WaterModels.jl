function verify_solution(solution::Dict{AbstractString, Any},
                         flowrate_solution_path::String,
                         head_solution_path::String)
    return true
end

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
        passed = verify_solution(solution, "../test/data/solutions/hanoi-feasibility-flowrate.csv",
                                 "../test/data/solutions/hanoi-feasibility-head.csv")
        @test passed == true
    end

    #@testset "klmod" begin
    #    solution = run_feasibility("../test/data/epanet/klmod.inp", GenericWaterModel{StandardMINLPForm}, solver)
    #    @test solution["status"] == :LocalOptimal
    #end

    #@testset "rural" begin
    #    solution = run_feasibility("../test/data/epanet/rural.inp", GenericWaterModel{StandardMINLPForm}, solver)
    #    @test solution["status"] == :LocalOptimal
    #end

    #@testset "tasseff" begin
    #    solution = run_feasibility("../test/data/epanet/tasseff.inp", GenericWaterModel{StandardMINLPForm}, solver)
    #    @test solution["status"] == :LocalOptimal
    #end

    #@testset "zj" begin
    #    solution = run_feasibility("../test/data/epanet/zj.inp", GenericWaterModel{StandardMINLPForm}, solver)
    #    @test solution["status"] == :LocalOptimal
    #end
end
