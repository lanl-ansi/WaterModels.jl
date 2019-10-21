@testset "Warm Start Functionality" begin
    @testset "Shamir network, cwf, CNLP formulation." begin
        network_data = WaterModels.parse_file("../test/data/epanet/shamir.inp")
        cold_result = WaterModels.run_cwf(network_data, CNLPWaterModel, ipopt)
        @test cold_result["termination_status"] == LOCALLY_SOLVED
        InfrastructureModels.update_data!(network_data, cold_result["solution"])

        set_start_all!(network_data)
        result = WaterModels.run_cwf(network_data, CNLPWaterModel, ipopt_ws)
        @test result["termination_status"] == LOCALLY_SOLVED
    end

    @testset "Shamir network, cwf, NCNLP formulation." begin
        network_data = WaterModels.parse_file("../test/data/epanet/shamir.inp")
        cold_result = WaterModels.run_cwf(network_data, NCNLPWaterModel, ipopt)
        @test cold_result["termination_status"] == LOCALLY_SOLVED
        InfrastructureModels.update_data!(network_data, cold_result["solution"])

        set_start_all!(network_data)
        result = WaterModels.run_cwf(network_data, NCNLPWaterModel, ipopt_ws)
        @test result["termination_status"] == LOCALLY_SOLVED
    end
end
