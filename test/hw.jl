@testset "Hazen-Williams MINLP Problems" begin
    @testset "Shamir network (physical feasibility, relaxation)." begin
        network_path = "../test/data/epanet/shamir.inp"
        solution = run_wf_hw(network_path, MINLPWaterModel, pavito)
        status = solution["status"]
        @test status == :LocalOptimal || status == :Optimal
    end

    @testset "Shamir network (diameter selection (reduced), relaxation)." begin
        network_path = "../test/data/epanet/shamir.inp"
        modification_path = "../test/data/json/shamir-reduced.json"
        solution = run_ne_hw(network_path, modification_path, MINLPWaterModel, pavito)
        status = solution["status"]
        @test status == :LocalOptimal || status == :Optimal
    end

    @testset "Shamir network (diameter selection (reduced), global algorithm)." begin
        network_path = "../test/data/epanet/shamir.inp"
        modification_path = "../test/data/json/shamir-reduced.json"
        status = solve_global(network_path, modification_path, ipopt, glpk)
        @test status == :LocalOptimal || status == :Optimal
    end
end
