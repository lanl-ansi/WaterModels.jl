@testset "Hazen-Williams CVXNLP Water Flow Problems" begin
    network_path = "../test/data/epanet/shamir.inp"

    @testset "Shamir network (unknown flow directions), CVXNLP formulation." begin
        solution = run_wf(network_path, CVXNLPWaterModel, ipopt, alpha=0.852)
        @test solution["termination_status"] == MOI.LOCALLY_SOLVED
    end

    @testset "Shamir network (unknown flow directions), MILP-R formulation." begin
        solution = run_wf(network_path, MILPRWaterModel, cbc, alpha=0.852)
        @test solution["termination_status"] == MOI.OPTIMAL
    end
end
