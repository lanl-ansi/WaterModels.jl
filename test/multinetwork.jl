@testset "replicate" begin
    shamir_mn_data = build_mn_data("../test/data/epanet/shamir.inp")

    @testset "Replicated Shamir network, WF problem, CNLP formulation." begin

        wm = build_generic_model(shamir_mn_data, CNLPWaterModel, WaterModels.post_wf, multinetwork=true)
        solution = solve_generic_model(wm, ipopt)
    end
end
