@testset "Warm Start Functionality" begin
    @testset "Shamir network, wf, NC formulation." begin
        data = WaterModels.parse_file("../examples/data/epanet/shamir.inp")
        cold_result = WaterModels.solve_wf(data, NCWaterModel, nlp_solver)
        @test _is_valid_status(cold_result["termination_status"])

        _IM.update_data!(data, cold_result["solution"])
        WaterModels.set_start_all!(data)
        WaterModels.fix_all_indicators!(data)
        WaterModels.fix_all_flow_directions!(data)

        result = WaterModels.solve_wf(data, NCWaterModel, nlp_solver)
        @test _is_valid_status(result["termination_status"])
    end
end
