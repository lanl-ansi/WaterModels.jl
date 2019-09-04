@testset "replicate" begin
    network_path_shamir_ts = "../test/data/epanet/shamir-ts.inp"

    @testset "Replicated Shamir network, WF problem, CNLP formulation." begin
        network_data = WaterModels.parse_file(network_path_shamir_ts)
        mn_data = WaterModels.make_multinetwork(network_data)
        wm = build_model(mn_data, CNLPWaterModel, WaterModels.post_mn_wf, multinetwork=true)
        solution = optimize_model!(wm, ipopt)
    end
end
