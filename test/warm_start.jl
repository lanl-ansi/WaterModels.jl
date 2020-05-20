@testset "Warm Start Functionality" begin
    @testset "Shamir network, wf, NLP formulation." begin
        data = WaterModels.parse_file("../test/data/epanet/shamir.inp")
        cold_result = WaterModels.solve_wf(data, NLPWaterModel, ipopt)
        @test cold_result["termination_status"] == LOCALLY_SOLVED
        InfrastructureModels.update_data!(data, cold_result["solution"])

        set_start_all!(data)
        result = WaterModels.solve_wf(data, NLPWaterModel, ipopt_ws)
        @test result["termination_status"] == LOCALLY_SOLVED
    end
end
