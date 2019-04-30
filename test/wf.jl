@testset "Hazen-Williams Water Flow Problems" begin
    network_path = "../test/data/epanet/shamir.inp"

    @testset "Shamir network (unknown flow directions), CNLP formulation." begin
        solution = run_wf(network_path, CNLPWaterModel, ipopt, alpha=1.852)
        @test solution["termination_status"] == MOI.LOCALLY_SOLVED
    end

    @testset "Shamir network (unknown flow directions), NCNLP formulation." begin
        solution = run_wf(network_path, NCNLPWaterModel, ipopt, alpha=1.852)
        @test solution["termination_status"] == MOI.LOCALLY_SOLVED
    end

    #@testset "Shamir network (unknown flow directions), MICP formulation." begin
    #    solution = run_wf(network_path, MICPWaterModel, juniper, alpha=1.852)
    #    @test solution["termination_status"] == MOI.LOCALLY_SOLVED
    #end
end
