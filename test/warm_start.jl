@testset "Warm Start Functionality" begin
    # Set Hazen-Williams-compatible network path.
    hw_network_path = "../test/data/epanet/shamir.inp"

    @testset "Shamir network (unknown flow directions), CNLP formulation." begin
        # This should take some number of iterations to converge.
        hw_network = WaterModels.parse_file(hw_network_path)
        cold_result = WaterModels.run_wf(hw_network, CNLPWaterModel, ipopt, alpha=1.852)
        InfrastructureModels.update_data!(hw_network, cold_result["solution"])

        # This should take fewer iterations to converge.
        WaterModels.set_start_directed_flow_rate!(hw_network)
        warm_result = WaterModels.run_wf(hw_network, CNLPWaterModel, ipopt_ws, alpha=1.852)
    end

    @testset "Shamir network (unknown flow directions), NCNLP formulation." begin
        # This should take some number of iterations to converge.
        hw_network = WaterModels.parse_file(hw_network_path)
        cold_result = WaterModels.run_wf(hw_network, NCNLPWaterModel, ipopt, alpha=1.852)
        InfrastructureModels.update_data!(hw_network, cold_result["solution"])

        # This should take fewer iterations to converge.
        WaterModels.set_start_head!(hw_network)
        WaterModels.set_start_undirected_flow_rate!(hw_network)
        warm_result = WaterModels.run_wf(hw_network, NCNLPWaterModel, ipopt_ws, alpha=1.852)
    end
end
