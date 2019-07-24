@testset "Optimal Water Flow Problems" begin
    # Set test network paths.
    richmond_skeleton_path = "../test/data/epanet/richmond-skeleton.inp"
    richmond_skeleton_sp_path = "../test/data/epanet/richmond-skeleton-sp.inp"

    @testset "Richmond network, single time point, NCNLP formulation." begin
        solution = run_owf(richmond_skeleton_sp_path, NCNLPWaterModel, ipopt, relaxed=true)
        #@test solution["termination_status"] == MOI.LOCALLY_SOLVED
    end

    @testset "Richmond network, multinetwork, NCNLP formulation." begin
        solution = run_mn_owf(richmond_skeleton_path, NCNLPWaterModel, ipopt, relaxed=true)
        #@test solution["termination_status"] == MOI.LOCALLY_SOLVED
    end
end
