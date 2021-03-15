@testset "src/core/data.jl" begin
    @testset "make_per_unit!" begin
        data = WaterModels.parse_file("../examples/data/json/shamir.json")
        make_per_unit!(data) # This function is currently experimental.
        @test data["per_unit"] == true
    end

    @testset "make_temporally_aggregated_multinetwork" begin
        network = WaterModels.parse_file("../test/data/epanet/multinetwork/owf-hw-lps.inp")
        network_mn = WaterModels.make_multinetwork(network)
        network_agg = make_temporally_aggregated_multinetwork(network_mn, [["1", "2"], ["3"]])

        @test network_mn["duration"] == network_agg["duration"]
        demand_old = network_mn["nw"]["1"]["demand"]["2"]["flow_nominal"] 
        @test demand_old < network_agg["nw"]["1"]["demand"]["2"]["flow_nominal"]

        network = WaterModels.parse_file("../test/data/epanet/multinetwork/prv-hw-lps.inp")
        network_mn = WaterModels.make_multinetwork(network)
        network_agg = make_temporally_aggregated_multinetwork(network_mn, [["1", "2"], ["3"]])

        @test network_mn["duration"] == network_agg["duration"]
        @test length(network_agg["nw"]) == 2
    end

    @testset "set_start! (single network)" begin
        data = WaterModels.parse_file("../test/data/epanet/snapshot/pipe-hw-lps.inp")
        data["pipe"]["1"]["q"] = 1.0 # Set the flow along the pipe.
        WaterModels.set_start!(data, "pipe", "q", "q_pipe_start")
    end

    @testset "set_start! (multinetwork)" begin
        data = WaterModels.parse_file("../test/data/epanet/snapshot/pipe-hw-lps.inp")
        data["pipe"]["1"]["q"] = 1.0 # Set the flow along the pipe.
        mn_data = WaterModels.replicate(data, 3)
        WaterModels.set_start!(mn_data, "pipe", "q", "q_pipe_start")
    end

    @testset "replicate" begin
        data = WaterModels.parse_file("../test/data/epanet/snapshot/pipe-hw-lps.inp")
        mn_data = WaterModels.replicate(data, 3)
        @test length(mn_data["nw"]) == 3
    end

    @testset "make_multinetwork shamir" begin
        network_data = WaterModels.parse_file("../test/data/epanet/multinetwork/shamir-ts.inp")

        @test !_IM.ismultinetwork(network_data)
        @test haskey(network_data, "time_series")
        @test isapprox(network_data["time_series"]["demand"]["2"]["flow_nominal"][1], 0.02777, rtol=1.0e-4)
        @test isapprox(network_data["time_series"]["demand"]["2"]["flow_nominal"][2], 0.5*0.02777, rtol=1.0e-4)
        @test isapprox(network_data["time_series"]["demand"]["2"]["flow_nominal"][3], 0.25*0.02777, rtol=1.0e-4)

        ts_length = network_data["time_series"]["num_steps"]
        mn_data = WaterModels.make_multinetwork(network_data)

        @test _IM.ismultinetwork(mn_data)
        @test !haskey(mn_data, "time_series")
        @test length(mn_data["nw"]) == ts_length
        @test isapprox(mn_data["nw"]["1"]["demand"]["2"]["flow_nominal"], 0.02777, rtol=1.0e-4)
        @test isapprox(mn_data["nw"]["2"]["demand"]["2"]["flow_nominal"], 0.5*0.02777, rtol=1.0e-4)
        @test isapprox(mn_data["nw"]["3"]["demand"]["2"]["flow_nominal"], 0.25*0.02777, rtol=1.0e-4)
    end

    @testset "split_multinetwork" begin
        data = WaterModels.parse_file("../test/data/epanet/snapshot/pipe-hw-lps.inp")
        mn_data = WaterModels.replicate(data, 3)
        new_mns = WaterModels.split_multinetwork(mn_data, [["1", "2"], ["2", "3"]])

        @test sort(collect(keys(new_mns[1]["nw"]))) == ["1", "2"]
        @test sort(collect(keys(new_mns[2]["nw"]))) == ["2", "3"]
    end

    @testset "epanet_to_watermodels!" begin
        network_data = WaterModels.parse_epanet("../test/data/epanet/snapshot/short-pipe-lps.inp")
        WaterModels.epanet_to_watermodels!(network_data)
        @test haskey(network_data["short_pipe"], "1")

        network_data = WaterModels.parse_epanet("../test/data/epanet/snapshot/short-pipe-valve-lps.inp")
        WaterModels.epanet_to_watermodels!(network_data)
        @test haskey(network_data["valve"], "1")
    end

    @testset "_read_status!" begin
        network_data = WaterModels.parse_epanet("../test/data/epanet/snapshot/shutoff_valve-status-hw-lps.inp")
        @test network_data["pipe"]["1"]["has_valve"]
        WaterModels.epanet_to_watermodels!(network_data)
        @test haskey(network_data["valve"], "1")
    end

    @testset "_read_controls!" begin
        network_data = WaterModels.parse_epanet("../test/data/epanet/snapshot/shutoff_valve-controls-hw-lps.inp")
        @test network_data["pipe"]["1"]["has_valve"]
        WaterModels.epanet_to_watermodels!(network_data)
        @test haskey(network_data["valve"], "1")
    end

    @testset "_remove_last_networks! (multinetwork)" begin
        data = WaterModels.parse_file("../test/data/epanet/snapshot/pipe-hw-lps.inp")
        mn_data = WaterModels.replicate(data, 3)
        WaterModels._remove_last_networks!(mn_data; last_num_steps = 1)
        @test length(mn_data["nw"]) == 2
    end
end
