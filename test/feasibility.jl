using CSV

function verify(solution::Dict{AbstractString, Any}, flowrate_solution_path::String, head_solution_path::String)
    # Parse the known head solution and compute the relative error.
    head = CSV.read(head_solution_path)
    nodes = merge(solution["solution"]["junctions"], solution["solution"]["reservoirs"])
    head[:test] = map(id -> nodes[string(id)]["h"], head[:id])
    head[:abs_err] = abs.(head[:solution] - head[:test])
    head_max_abs_err = maximum(head[:abs_err])

    # Parse the known flow solution and compute the relative error.
    flow = CSV.read(flowrate_solution_path)
    edges = solution["solution"]["pipes"]
    flow[:test] = map(id -> edges[string(id)]["q"], flow[:id])
    flow[:abs_err] = abs.(flow[:solution] - flow[:test])
    flow_max_abs_err = maximum(flow[:abs_err])

    # If the maximum absolute error is small, then the solution is valid.
    if head_max_abs_err < 0.10 && flow_max_abs_err < 0.10
        return true
    else
        return false
    end
end

@testset "Hazen-Williams MINLP Problems" begin
    @testset "foss_poly_1" begin
        solution = run_feasibility("../test/data/epanet/foss_poly_1.inp", GenericWaterModel{StandardMINLPForm}, solver)
        @test solution["status"] == :LocalOptimal
        # @test verify(solution, "../test/data/solutions/foss_poly_1-feasibility-flowrate.csv", "../test/data/solutions/foss_poly_1-feasibility-head.csv")
    end

    @testset "hanoi_extended" begin
        solution = run_feasibility("../test/data/epanet/hanoi_extended.inp", GenericWaterModel{StandardMINLPForm}, solver)
        @test solution["status"] == :LocalOptimal
        # @test verify(solution, "../test/data/solutions/hanoi_extended-feasibility-flowrate.csv", "../test/data/solutions/hanoi_extended-feasibility-head.csv")
    end

    @testset "hanoi" begin
        solution = run_feasibility("../test/data/epanet/hanoi.inp", GenericWaterModel{StandardMINLPForm}, solver)
        @test solution["status"] == :LocalOptimal
        # @test verify(solution, "../test/data/solutions/hanoi-feasibility-flowrate.csv", "../test/data/solutions/hanoi-feasibility-head.csv")
    end

    @testset "zj" begin
        solution = run_feasibility("../test/data/epanet/zj.inp", GenericWaterModel{StandardMINLPForm}, solver)
        @test solution["status"] == :LocalOptimal
        # @test solution["status"] == :LocalInfeasible # The actual EPANET solution is infeasible.
    end

    # This one is too large to solve in a unit test.
    # @testset "klmod" begin
    #     solution = run_feasibility("../test/data/epanet/klmod.inp", GenericWaterModel{StandardMINLPForm}, solver)
    #     @test solution["status"] == :LocalOptimal
    #     # @test verify(solution, "../test/data/solutions/klmod-feasibility-flowrate.csv", "../test/data/solutions/klmod-feasibility-head.csv")
    # end
end

@testset "Darcy-Weisbach MINLP Problems" begin
    @testset "balerma" begin
        solution = run_feasibility("../test/data/epanet/balerma.inp", GenericWaterModel{StandardMINLPForm}, solver)
        @test solution["status"] == :LocalOptimal
        # @test verify(solution, "../test/data/solutions/balerma-feasibility-flowrate.csv", "../test/data/solutions/balerma-feasibility-head.csv")
    end

    # TODO: Why doesn't a feasible solution exist to this problem?
    # @testset "rural" begin
    #     solution = run_feasibility("../test/data/epanet/rural.inp", GenericWaterModel{StandardMINLPForm}, solver)
    #     @test solution["status"] == :LocalOptimal
    #     # @test verify(solution, "../test/data/solutions/rural-feasibility-flowrate.csv", "../test/data/solutions/rural-feasibility-head.csv")
    # end
end
