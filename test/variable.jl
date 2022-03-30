@testset "src/core/variable.jl" begin
    @testset "variable_reservoir_flow" begin
        network = WaterModels.parse_file("../test/data/epanet/snapshot/pipe-hw-lps.inp")
        pipe = network["pipe"]["1"] # Get pipe that will be reversed.
        pipe["node_fr"], pipe["node_to"] = pipe["node_to"], pipe["node_fr"]
        wm = WaterModels.instantiate_model(network, WaterModels.LRDWaterModel, _build_null_model)
        WaterModels.variable_reservoir_flow(wm)
        @test isa(WaterModels.var(wm, :q_reservoir, 1), JuMP.VariableRef)
    end

    @testset "variable_valve_indicator (relax = true)" begin
        network = WaterModels.parse_file("../test/data/epanet/snapshot/shutoff_valve-hw-lps.inp")
        wm = WaterModels.instantiate_model(network, WaterModels.LRDWaterModel, _build_null_model)
        WaterModels.variable_valve_indicator(wm; relax = true)
        @test !JuMP.is_binary(WaterModels.var(wm, :z_valve, 1))
        @test JuMP.lower_bound(WaterModels.var(wm, :z_valve, 1)) == 0.0
        @test JuMP.upper_bound(WaterModels.var(wm, :z_valve, 1)) == 1.0
    end

    @testset "variable_regulator_indicator (relax = true)" begin
        network = WaterModels.parse_file("../test/data/epanet/snapshot/prv-hw-lps.inp")
        wm = WaterModels.instantiate_model(network, WaterModels.LRDWaterModel, _build_null_model)
        WaterModels.variable_regulator_indicator(wm; relax = true)
        @test !JuMP.is_binary(WaterModels.var(wm, :z_regulator, 2))
        @test JuMP.lower_bound(WaterModels.var(wm, :z_regulator, 2)) == 0.0
        @test JuMP.upper_bound(WaterModels.var(wm, :z_regulator, 2)) == 1.0
    end

    @testset "variable_pump_indicator (relax = true)" begin
        network = WaterModels.parse_file("../test/data/epanet/snapshot/pump-hw-lps.inp")
        wm = WaterModels.instantiate_model(network, WaterModels.LRDWaterModel, _build_null_model)
        WaterModels.variable_pump_indicator(wm; relax = true)
        @test !JuMP.is_binary(WaterModels.var(wm, :z_pump, 1))
        @test JuMP.lower_bound(WaterModels.var(wm, :z_pump, 1)) == 0.0
        @test JuMP.upper_bound(WaterModels.var(wm, :z_pump, 1)) == 1.0
    end

    @testset "variable_pump_switch_on (relax = true)" begin
        network = WaterModels.parse_file("../test/data/epanet/snapshot/pump-hw-lps.inp")
        wm = WaterModels.instantiate_model(network, WaterModels.LRDWaterModel, _build_null_model)
        WaterModels.variable_pump_switch_on(wm; relax = true)
        @test !JuMP.is_binary(WaterModels.var(wm, :z_switch_on_pump, 1))
        @test JuMP.lower_bound(WaterModels.var(wm, :z_switch_on_pump, 1)) == 0.0
        @test JuMP.upper_bound(WaterModels.var(wm, :z_switch_on_pump, 1)) == 1.0
    end

    @testset "variable_pump_switch_off (relax = true)" begin
        network = WaterModels.parse_file("../test/data/epanet/snapshot/pump-hw-lps.inp")
        wm = WaterModels.instantiate_model(network, WaterModels.LRDWaterModel, _build_null_model)
        WaterModels.variable_pump_switch_off(wm; relax = true)
        @test !JuMP.is_binary(WaterModels.var(wm, :z_switch_off_pump, 1))
        @test JuMP.lower_bound(WaterModels.var(wm, :z_switch_off_pump, 1)) == 0.0
        @test JuMP.upper_bound(WaterModels.var(wm, :z_switch_off_pump, 1)) == 1.0
    end

    @testset "variable_des_pipe_indicator (relax = true)" begin
        network = parse_file("../test/data/json/shamir.json")
        wm = WaterModels.instantiate_model(network, WaterModels.LRDWaterModel, _build_null_model)
        WaterModels.variable_des_pipe_indicator(wm; relax = true)
        @test !JuMP.is_binary(WaterModels.var(wm, :z_des_pipe, 1))
        @test JuMP.lower_bound(WaterModels.var(wm, :z_des_pipe, 1)) == 0.0
        @test JuMP.upper_bound(WaterModels.var(wm, :z_des_pipe, 1)) == 1.0
    end

    @testset "_fix_indicator_variable" begin
        @testset "fix variable to one" begin
            network = WaterModels.parse_file("../test/data/epanet/snapshot/pump-hw-lps.inp")
            network["pump"]["1"]["z_min"], network["pump"]["1"]["z_max"] = 1.0, 1.0
            wm = WaterModels.instantiate_model(network, WaterModels.LRDWaterModel, _build_null_model)
            WaterModels.variable_pump_indicator(wm; relax = false)
            WaterModels._fix_indicator_variable(WaterModels.var(wm, :z_pump, 1), network["pump"]["1"], "z")
            @test JuMP.is_fixed(WaterModels.var(wm, :z_pump, 1))
            @test JuMP.fix_value(WaterModels.var(wm, :z_pump, 1)) == 1.0
        end

        @testset "fix variable to zero" begin
            network = WaterModels.parse_file("../test/data/epanet/snapshot/pump-hw-lps.inp")
            network["pump"]["1"]["z_min"], network["pump"]["1"]["z_max"] = 0.0, 0.0
            wm = WaterModels.instantiate_model(network, WaterModels.LRDWaterModel, _build_null_model)
            WaterModels.variable_pump_indicator(wm; relax = false)
            WaterModels._fix_indicator_variable(WaterModels.var(wm, :z_pump, 1), network["pump"]["1"], "z")
            @test JuMP.is_fixed(WaterModels.var(wm, :z_pump, 1))
            @test JuMP.fix_value(WaterModels.var(wm, :z_pump, 1)) == 0.0
        end
    end
end