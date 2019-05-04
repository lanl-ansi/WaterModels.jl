@testset "Hazen-Williams MICP Problems" begin
    @testset "Shamir network (physical feasibility, MICP)." begin
        network_path = "../test/data/epanet/shamir.inp"
        result = run_wf_hw(network_path, MICPWaterModel, pavito)
        @test result["status"] == :LocalOptimal || result["status"] == :Optimal
    end

    @testset "Shamir network (diameter selection (reduced), MICP)." begin
        network_path = "../test/data/epanet/shamir.inp"
        modification_path = "../test/data/json/shamir-reduced.json"
        result = run_ne_hw(network_path, modification_path, MICPWaterModel, pavito)
        @test result["status"] == :LocalOptimal || result["status"] == :Optimal
    end

    @testset "Shamir network (diameter selection (reduced), global algorithm)." begin
        network_path = "../test/data/epanet/shamir.inp"
        modification_path = "../test/data/json/shamir.json"
        status = solve_global(network_path, modification_path, ipopt, glpk)
        @test status == :LocalOptimal || status == :Optimal
    end
end
