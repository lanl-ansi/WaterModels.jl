@testset "Warm Start Functionality" begin
    @testset "Shamir network, wf, NC formulation." begin
        data = WaterModels.parse_file("../examples/data/epanet/shamir.inp")
        cold_result = WaterModels.run_wf(data, NCWaterModel, ipopt)
        @test _is_valid_status(cold_result["termination_status"])

        _IM.update_data!(data, cold_result["solution"])
        WaterModels.set_start_all!(data)
        WaterModels.fix_all_indicators!(data)

        result = WaterModels.run_wf(data, NCWaterModel, ipopt)
        @test _is_valid_status(result["termination_status"])
    end
end
