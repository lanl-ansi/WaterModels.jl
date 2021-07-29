@testset "src/util/obbt.jl" begin
    @testset "solve_obbt_owf! (pipe)" begin
        data = WaterModels.parse_file("../test/data/epanet/multinetwork/pipe-hw-lps.inp")
        set_flow_partitions!(data, 1.0, 1.0e-4)
        
        solve_obbt_owf!(data, cbc; model_type = LRDWaterModel, max_iter = 1)
        @test haskey(data["pipe"]["1"], "flow_min")
        @test haskey(data["pipe"]["1"], "flow_max")
    end

    @testset "solve_obbt_owf! (pump)" begin
        data = WaterModels.parse_file("../test/data/epanet/multinetwork/pump-hw-lps.inp")
        set_flow_partitions!(data, 1.0, 1.0e-4)

        solve_obbt_owf!(data, cbc; model_type = LRDWaterModel, max_iter = 1)
        @test haskey(data["pump"]["1"], "flow_max")
    end

    @testset "solve_obbt_owf! (regulator)" begin
        data = WaterModels.parse_file("../test/data/epanet/multinetwork/prv-hw-lps.inp")
        set_flow_partitions!(data, 1.0, 1.0e-4)
        solve_obbt_owf!(data, cbc; model_type = LRDWaterModel, max_iter = 1)
        @test haskey(data["regulator"]["2"], "flow_max")
    end

    @testset "solve_obbt_owf! (short pipe)" begin
        data = WaterModels.parse_file("../test/data/epanet/multinetwork/short-pipe-lps.inp")
        WaterModels.convert_short_pipes!(data)
        set_flow_partitions!(data, 1.0, 1.0e-4)

        solve_obbt_owf!(data, cbc; model_type = LRDWaterModel, max_iter = 1)
        @test haskey(data["short_pipe"]["1"], "flow_min")
        @test haskey(data["short_pipe"]["1"], "flow_max")
    end

    @testset "solve_obbt_owf! (valve)" begin
        data = WaterModels.parse_file("../test/data/epanet/multinetwork/shutoff_valve-hw-lps.inp")
        set_flow_partitions!(data, 1.0, 1.0e-4)

        solve_obbt_owf!(data, cbc; model_type = LRDWaterModel, max_iter = 1)
        @test haskey(data["valve"]["1"], "flow_min")
        @test haskey(data["valve"]["1"], "flow_max")
    end
end
