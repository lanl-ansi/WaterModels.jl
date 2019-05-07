@testset "Warm Start Functionality" begin
    # Set Hazen-Williams-compatible network path.
    hw_network_path = "../test/data/epanet/shamir.inp"

    @testset "Shamir network, water flow, CNLP formulation." begin
        # This should take some number of iterations to converge.
        hw_network = WaterModels.parse_file(hw_network_path)
        cold_result = WaterModels.run_wf(hw_network, CNLPWaterModel, ipopt, alpha=1.852)
        @test cold_result["termination_status"] == MOI.LOCALLY_SOLVED
        InfrastructureModels.update_data!(hw_network, cold_result["solution"])

        # This should take fewer iterations to converge.
        WaterModels.set_start_directed_flow_rate!(hw_network)
        result = WaterModels.run_wf(hw_network, CNLPWaterModel, ipopt_ws, alpha=1.852)
        @test result["termination_status"] == MOI.LOCALLY_SOLVED
    end

    @testset "Shamir network, water flow, NCNLP formulation." begin
        # This should take some number of iterations to converge.
        hw_network = WaterModels.parse_file(hw_network_path)
        cold_result = WaterModels.run_wf(hw_network, NCNLPWaterModel, ipopt, alpha=1.852)
        @test cold_result["termination_status"] == MOI.LOCALLY_SOLVED
        InfrastructureModels.update_data!(hw_network, cold_result["solution"])

        # This should take fewer iterations to converge.
        WaterModels.set_start_head!(hw_network)
        WaterModels.set_start_undirected_flow_rate!(hw_network)
        result = WaterModels.run_wf(hw_network, NCNLPWaterModel, ipopt_ws, alpha=1.852)
        @test result["termination_status"] == MOI.LOCALLY_SOLVED
    end

    @testset "Shamir network, expansion planning, NCNLP formulation." begin
        # This should take some number of iterations to converge.
        hw_network = WaterModels.parse_file(hw_network_path)
        cold_result = WaterModels.run_wf(hw_network, NCNLPWaterModel, ipopt, alpha=1.852)
        @test cold_result["termination_status"] == MOI.LOCALLY_SOLVED
        InfrastructureModels.update_data!(hw_network, cold_result["solution"])

        # This should take fewer iterations to converge.
        modifications = WaterModels.parse_file("../test/data/json/shamir.json")
        InfrastructureModels.update_data!(hw_network, modifications)
        WaterModels.set_start_head!(hw_network)
        WaterModels.set_start_undirected_flow_rate!(hw_network)
        WaterModels.set_start_undirected_flow_rate_ne!(hw_network)
        WaterModels.set_start_resistance_ne!(hw_network)
        result = WaterModels.run_ne(hw_network, NCNLPWaterModel, ipopt_ws, alpha=1.852, relaxed=true)
        @test result["termination_status"] == MOI.LOCALLY_SOLVED
    end
end
