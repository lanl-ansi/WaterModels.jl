@testset "src/util/variable_index.jl" begin
    data = WaterModels.parse_file("../test/data/epanet/multinetwork/owf-hw-lps.inp")
    mn_data = WaterModels.make_multinetwork(data)
    wm = instantiate_model(mn_data, CRDWaterModel, build_mn_wf)

    @testset "_VariableIndex instantiation" begin
        vid = WaterModels._VariableIndex(1, :pump, :z_pump, 1)
        @test vid.network_index == 1
        @test vid.component_type == :pump
        @test vid.variable_symbol == :z_pump
        @test vid.component_index == 1
    end

    @testset "_get_variable_from_index" begin
        vid = WaterModels._VariableIndex(1, :pump, :z_pump, 1)
        z_pump = WaterModels._get_variable_from_index(wm, vid)
        @test JuMP.is_binary(z_pump) == true
    end

    @testset "_get_lower_bound_from_index (JuMP.VariableRef)" begin
        vid = WaterModels._VariableIndex(1, :pump, :z_pump, 1)
        @test WaterModels._get_lower_bound_from_index(wm, vid) == 0.0
    end

    @testset "_get_lower_bound_from_index (JuMP.AffExpr)" begin
        vid = WaterModels._VariableIndex(1, :pipe, :q_pipe, 2)
        @test WaterModels._get_lower_bound_from_index(wm, vid) < 0.0
    end

    @testset "_get_upper_bound_from_index (JuMP.VariableRef)" begin
        vid = WaterModels._VariableIndex(1, :pump, :z_pump, 1)
        @test WaterModels._get_upper_bound_from_index(wm, vid) == 1.0
    end

    @testset "_get_upper_bound_from_index (JuMP.AffExpr)" begin
        vid = WaterModels._VariableIndex(1, :pipe, :q_pipe, 2)
        @test WaterModels._get_upper_bound_from_index(wm, vid) > 0.0
    end

    @testset "_get_indicator_variable_indices" begin
        vids = WaterModels._get_indicator_variable_indices(wm; nw = 1)
        vid = WaterModels._VariableIndex(1, :pump, :z_pump, 1)
        @test !any(x -> x.variable_symbol == :y_pump, vids)
        @test any(x -> x.variable_symbol == :z_pump, vids)
    end

    @testset "_get_direction_variable_indices (AbstractDirectedModel)" begin
        vids = WaterModels._get_direction_variable_indices(wm; nw=1)
        @test any(x -> x.variable_symbol == :y_pump, vids)
        @test !any(x -> x.variable_symbol == :z_pump, vids)
    end

    @testset "_get_binary_variable_indices (AbstractDirectedModel)" begin
        vids = WaterModels._get_binary_variable_indices(wm; nw=1)
        @test any(x -> x.variable_symbol == :y_pump, vids)
        @test any(x -> x.variable_symbol == :z_pump, vids)
    end

    @testset "_get_binary_variable_indices (AbstractUndirectedModel)" begin
        wm_ud = instantiate_model(mn_data, NCWaterModel, build_mn_wf)
        vids = WaterModels._get_binary_variable_indices(wm_ud; nw=1)
        @test !any(x -> x.variable_symbol == :y_pump, vids)
        @test any(x -> x.variable_symbol == :z_pump, vids)
    end
end
