@testset "src/core/data.jl" begin
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
        network_data = WaterModels.parse_file("../test/data/epanet/shamir-ts.inp")

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
end
