@testset "src/util/obbt.jl" begin
    @testset "solve_obbt_owf! (pipe)" begin
        data = WaterModels.parse_file("../test/data/epanet/multinetwork/pipe-hw-lps.inp")
        set_flow_partitions_si!(data, 10.0, 1.0e-4)
        solve_obbt_owf!(data, cbc; model_type = LRDWaterModel, max_iter = 2, solve_relaxed = true)
        @test haskey(data["pipe"]["1"], "flow_min")
        @test haskey(data["pipe"]["1"], "flow_max")

        network_mn = WaterModels.make_multinetwork(data)
        set_flow_partitions_si!(network_mn, 10.0, 1.0e-4)
        solve_obbt_owf!(network_mn, cbc; model_type = LRDWaterModel,
            max_iter = 2, use_relaxed_network = false, solve_relaxed = true)
        @test haskey(network_mn["nw"]["1"]["pipe"]["1"], "flow_min")
        @test haskey(network_mn["nw"]["1"]["pipe"]["1"], "flow_max")
    end

    @testset "solve_obbt_owf! (pump)" begin
        data = WaterModels.parse_file("../test/data/epanet/multinetwork/pump-hw-lps.inp")
        set_flow_partitions_si!(data, 10.0, 1.0e-4)
        solve_obbt_owf!(data, cbc; model_type = LRDWaterModel, max_iter = 2, solve_relaxed = true)
        @test haskey(data["pump"]["1"], "flow_max")

        network_mn = WaterModels.make_multinetwork(data)
        set_flow_partitions_si!(network_mn, 10.0, 1.0e-4)
        solve_obbt_owf!(network_mn, cbc; model_type = LRDWaterModel,
            max_iter = 2, use_relaxed_network = false, solve_relaxed = true)
        @test haskey(network_mn["nw"]["1"]["pump"]["1"], "flow_max")
    end

    @testset "solve_obbt_owf! (regulator)" begin
        data = WaterModels.parse_file("../test/data/epanet/multinetwork/prv-hw-lps.inp")
        set_flow_partitions_si!(data, 10.0, 1.0e-4)
        solve_obbt_owf!(data, cbc; model_type = LRDWaterModel, max_iter = 2, solve_relaxed = true)
        @test haskey(data["regulator"]["2"], "flow_max")

        network_mn = WaterModels.make_multinetwork(data)
        set_flow_partitions_si!(network_mn, 10.0, 1.0e-4)
        solve_obbt_owf!(network_mn, cbc; model_type = LRDWaterModel,
            max_iter = 2, use_relaxed_network = false, solve_relaxed = true)
        @test haskey(network_mn["nw"]["2"]["regulator"]["2"], "flow_max")
    end

    @testset "solve_obbt_owf! (short pipe)" begin
        data = WaterModels.parse_file("../test/data/epanet/multinetwork/short-pipe-lps.inp")
        WaterModels.convert_short_pipes!(data)
        set_flow_partitions_si!(data, 10.0, 1.0e-4)
        solve_obbt_owf!(data, cbc; model_type = LRDWaterModel, max_iter = 2, solve_relaxed = true)
        @test haskey(data["short_pipe"]["1"], "flow_min")
        @test haskey(data["short_pipe"]["1"], "flow_max")

        network_mn = WaterModels.make_multinetwork(data)
        set_flow_partitions_si!(network_mn, 10.0, 1.0e-4)
        solve_obbt_owf!(network_mn, cbc; model_type = LRDWaterModel,
            max_iter = 2, use_relaxed_network = false, solve_relaxed = true)
        @test haskey(network_mn["nw"]["1"]["short_pipe"]["1"], "flow_min")
        @test haskey(network_mn["nw"]["1"]["short_pipe"]["1"], "flow_max")
    end

    @testset "solve_obbt_owf! (valve)" begin
        data = WaterModels.parse_file("../test/data/epanet/multinetwork/shutoff_valve-hw-lps.inp")
        set_flow_partitions_si!(data, 10.0, 1.0e-4)
        solve_obbt_owf!(data, cbc; model_type = LRDWaterModel, max_iter = 2, solve_relaxed = true)
        @test haskey(data["valve"]["1"], "flow_min")
        @test haskey(data["valve"]["1"], "flow_max")

        network_mn = WaterModels.make_multinetwork(data)
        set_flow_partitions_si!(network_mn, 10.0, 1.0e-4)
        solve_obbt_owf!(network_mn, cbc; model_type = LRDWaterModel,
            max_iter = 2, use_relaxed_network = false, solve_relaxed = true)
        @test haskey(network_mn["nw"]["1"]["valve"]["1"], "flow_min")
        @test haskey(network_mn["nw"]["1"]["valve"]["1"], "flow_max")
    end
end
