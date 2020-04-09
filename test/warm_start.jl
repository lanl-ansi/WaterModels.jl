@testset "Warm Start Functionality" begin
    @testset "Shamir network, wf, CNLP formulation." begin
        data = WaterModels.parse_file("../test/data/epanet/shamir.inp")
        cold_result = WaterModels.solve_wf(data, CNLPWaterModel, ipopt,
            setting=Dict("output"=>Dict("duals"=>true)),
            solution_processors=[sol_data_model!])
        @test cold_result["termination_status"] == LOCALLY_SOLVED
        InfrastructureModels.update_data!(data, cold_result["solution"])

        set_start_all!(data)
        result = WaterModels.solve_wf(data, CNLPWaterModel, ipopt_ws)
        @test result["termination_status"] == LOCALLY_SOLVED
    end

    @testset "Shamir network, wf, NCNLP formulation." begin
        data = WaterModels.parse_file("../test/data/epanet/shamir.inp")
        cold_result = WaterModels.solve_wf(data, NCNLPWaterModel, ipopt)
        @test cold_result["termination_status"] == LOCALLY_SOLVED
        InfrastructureModels.update_data!(data, cold_result["solution"])

        set_start_all!(data)
        result = WaterModels.solve_wf(data, NCNLPWaterModel, ipopt_ws)
        @test result["termination_status"] == LOCALLY_SOLVED
    end
end
