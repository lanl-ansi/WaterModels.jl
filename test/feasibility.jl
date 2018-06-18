using CSV

function verify(solution::Dict{AbstractString, Any}, flowrate_solution_path::String, head_solution_path::String)
    # Parse the known head solution and compute the relative error.
    head = CSV.read(headution_path)
    nodes = merge(solution["solution"]["junctions"], solution["solution"]["reservoirs"])
    head[:test] = map(id -> nodes[string(id)]["h"], head[:id])
    head[:rel_err] = abs((head[:solution] - head_sol[:test]) ./ head_sol[:solution])
    head_max_rel_err = max(head[:rel_err])

    # Parse the known flow solution and compute the relative error.
    flow = CSV.read(flowrate_solution_path)
    edges = solution["solution"]["pipes"]
    flow[:test] = map(id -> edges[string(id)]["q"], flow[:id])
    flow[:rel_err] = abs((flow[:solution] - flow[:test]) ./ flow[:solution])
    flow_max_rel_err = max(flow[:rel_err])

    # If the maximum relative error is small, then the solution is good.
    if head_max_rel_err < 1.0e-4 && flow_max_rel_err < 1.0e-4
        return true
    else
        return false
    end
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

    @testset "hanoi_extended" begin
        solution = run_feasibility("../test/data/epanet/hanoi_extended.inp", GenericWaterModel{StandardMINLPForm}, solver)
        @test solution["status"] == :LocalOptimal
        @test verify(solution, "../test/data/solutions/hanoi-extended-feasibility-flowrate.csv", "../test/data/solutions/hanoi-feasibility-head.csv")
    end

    @testset "hanoi" begin
        solution = run_feasibility("../test/data/epanet/hanoi.inp", GenericWaterModel{StandardMINLPForm}, solver)
        @test solution["status"] == :LocalOptimal
        @test verify(solution, "../test/data/solutions/hanoi-feasibility-flowrate.csv", "../test/data/solutions/hanoi-feasibility-head.csv")
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
