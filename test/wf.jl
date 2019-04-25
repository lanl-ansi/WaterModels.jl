@testset "Hazen-Williams Water Flow Problems" begin
    network_path = "../test/data/epanet/shamir.inp"

    @testset "Shamir network (unknown flow directions), NCNLP formulation." begin
        solution = run_wf(network_path, NCNLPWaterModel, ipopt, alpha=1.852)
        @test solution["termination_status"] == MOI.LOCALLY_SOLVED
    end

    #@testset "Shamir network (unknown flow directions), CVXNLP formulation." begin
    #    solution = run_wf(network_path, CVXNLPWaterModel, ipopt, alpha=1.852)
    #    @test solution["termination_status"] == MOI.LOCALLY_SOLVED
    #end
end
