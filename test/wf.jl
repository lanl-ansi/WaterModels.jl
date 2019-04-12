@testset "Hazen-Williams CVXNLP Problems" begin
    @testset "Shamir network (unknown flow directions)." begin
        network_path = "../test/data/epanet/shamir.inp"
        solution = run_wf(network_path, CVXNLPWaterModel, ipopt, alpha=0.852)
        @test solution["termination_status"] == MOI.LOCALLY_SOLVED
    end
end
