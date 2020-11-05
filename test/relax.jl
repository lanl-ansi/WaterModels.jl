@testset "src/util/relax.jl" begin
    @testset "relax_all_binary_variables!" begin
        data = WaterModels.parse_file("../test/data/epanet/multinetwork/owf-hw-lps.inp")
        mn_data = WaterModels.make_multinetwork(data)
        wm = instantiate_model(mn_data, CRDWaterModel, build_mn_wf)
        WaterModels.relax_all_binary_variables!(wm)

        z_pump_1 = WaterModels.var(wm, 1, :z_pump, 1)
        @test JuMP.is_binary(z_pump_1) == false
        z_pump_3 = WaterModels.var(wm, 3, :z_pump, 1)
        @test JuMP.is_binary(z_pump_3) == false
        y_pump_1 = WaterModels.var(wm, 1, :y_pump, 1)
        @test JuMP.is_binary(y_pump_1) == false
        y_pump_3 = WaterModels.var(wm, 3, :y_pump, 1)
        @test JuMP.is_binary(y_pump_3) == false
    end

    @testset "_relax_binary_variable! (unfixed variable)" begin
        data = WaterModels.parse_file("../test/data/epanet/multinetwork/owf-hw-lps.inp")
        mn_data = WaterModels.make_multinetwork(data)
        wm = instantiate_model(mn_data, NCWaterModel, build_mn_wf)
        z_pump = WaterModels.var(wm, 1, :z_pump, 1)
        WaterModels._relax_binary_variable!(z_pump)
        @test JuMP.is_binary(z_pump) == false
    end

    @testset "_relax_binary_variable! (fixed variable)" begin
        data = WaterModels.parse_file("../test/data/epanet/multinetwork/owf-hw-lps.inp")
        mn_data = WaterModels.make_multinetwork(data)
        wm = instantiate_model(mn_data, NCWaterModel, build_mn_wf)
        z_pump = WaterModels.var(wm, 1, :z_pump, 1)
        WaterModels._relax_binary_variable!(z_pump)
        @test JuMP.is_binary(z_pump) == false
    end

    @testset "_relax_variables_with_symbol!" begin
        data = WaterModels.parse_file("../test/data/epanet/multinetwork/owf-hw-lps.inp")
        mn_data = WaterModels.make_multinetwork(data)
        wm = instantiate_model(mn_data, NCWaterModel, build_mn_wf)
        WaterModels._relax_variables_with_symbol!(wm, :z_pump)

        z_pump_1 = WaterModels.var(wm, 1, :z_pump, 1)
        @test JuMP.is_binary(z_pump_1) == false
        z_pump_3 = WaterModels.var(wm, 3, :z_pump, 1)
        @test JuMP.is_binary(z_pump_3) == false
    end

    @testset "_relax_all_direction_variables!" begin
        data = WaterModels.parse_file("../test/data/epanet/multinetwork/owf-hw-lps.inp")
        mn_data = WaterModels.make_multinetwork(data)
        wm = instantiate_model(mn_data, CRDWaterModel, build_mn_wf)
        WaterModels._relax_all_direction_variables!(wm)

        y_pump_1 = WaterModels.var(wm, 1, :y_pump, 1)
        @test JuMP.is_binary(y_pump_1) == false
        y_pump_3 = WaterModels.var(wm, 3, :y_pump, 1)
        @test JuMP.is_binary(y_pump_3) == false
    end

    @testset "_relax_all_indicator_variables!" begin
        data = WaterModels.parse_file("../test/data/epanet/multinetwork/owf-hw-lps.inp")
        mn_data = WaterModels.make_multinetwork(data)
        wm = instantiate_model(mn_data, NCWaterModel, build_mn_wf)
        WaterModels._relax_all_indicator_variables!(wm)

        z_pump_1 = WaterModels.var(wm, 1, :z_pump, 1)
        @test JuMP.is_binary(z_pump_1) == false
        z_pump_3 = WaterModels.var(wm, 3, :z_pump, 1)
        @test JuMP.is_binary(z_pump_3) == false
    end

    @testset "_relax_indicator_variables!" begin
        data = WaterModels.parse_file("../test/data/epanet/multinetwork/owf-hw-lps.inp")
        mn_data = WaterModels.make_multinetwork(data)
        wm = instantiate_model(mn_data, NCWaterModel, build_mn_wf)
        WaterModels._relax_last_indicator_variables!(wm; last_num_steps = 2)

        z_pump_1 = WaterModels.var(wm, 1, :z_pump, 1)
        @test JuMP.is_binary(z_pump_1) == true
        z_pump_2 = WaterModels.var(wm, 2, :z_pump, 1)
        @test JuMP.is_binary(z_pump_2) == false
        z_pump_3 = WaterModels.var(wm, 3, :z_pump, 1)
        @test JuMP.is_binary(z_pump_3) == false
    end
end
