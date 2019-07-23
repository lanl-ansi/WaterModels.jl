@testset "Optimal Water Flow Problems" begin
    # Set test network paths.
    richmond_skeleton_path = "../test/data/epanet/richmond-skeleton.inp"

    @testset "Richmond network (unknown flow directions), NCNLP formulation." begin
        solution = run_mn_owf(richmond_skeleton_path, NCNLPWaterModel, ipopt)
    end
end
