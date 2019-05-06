@testset "Water Flow Problems" begin
    # Set Hazen-Williams-compatible network path.
    hw_network_path = "../test/data/epanet/shamir.inp"

    # Set Darcy-Weisbach-compatible network path.
    dw_network_path = "../test/data/epanet/balerma.inp"

    @testset "Shamir network (unknown flow directions), CNLP formulation." begin
        solution = run_wf(hw_network_path, CNLPWaterModel, ipopt, alpha=1.852)
        @test solution["termination_status"] == MOI.LOCALLY_SOLVED
        @test isapprox(solution["solution"]["pipes"]["2"]["q"], 0.093565, rtol=1.0e-4)
        @test isapprox(solution["solution"]["pipes"]["6"]["q"], 0.055710, rtol=1.0e-4)
    end

    @testset "Shamir network (unknown flow directions), NCNLP formulation." begin
        solution = run_wf(hw_network_path, NCNLPWaterModel, ipopt, alpha=1.852)
        @test solution["termination_status"] == MOI.LOCALLY_SOLVED
        @test isapprox(solution["solution"]["pipes"]["2"]["q"], 0.093565, rtol=1.0e-4)
        @test isapprox(solution["solution"]["pipes"]["6"]["q"], 0.055710, rtol=1.0e-4)
        @test isapprox(solution["solution"]["junctions"]["2"]["h"], 203.247650, rtol=1.0e-4)
        @test isapprox(solution["solution"]["junctions"]["6"]["h"], 195.445953, rtol=1.0e-4)
    end

    @testset "Shamir network (unknown flow directions), MICP formulation." begin
        solution = run_wf(hw_network_path, MICPWaterModel, ipopt, alpha=1.852, relaxed=true)
        @test solution["termination_status"] == MOI.LOCALLY_SOLVED
    end

    @testset "Balerma network (unknown flow directions), CNLP formulation." begin
        solution = run_wf(dw_network_path, CNLPWaterModel, ipopt, alpha=2.0)
        @test solution["termination_status"] == MOI.LOCALLY_SOLVED
    end
end
