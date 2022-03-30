@testset "src/core/data.jl" begin
    @testset "make_per_unit!" begin
        data = WaterModels.parse_file("../examples/data/json/shamir.json")
        make_per_unit!(data) # This function is currently experimental.
        @test data["per_unit"] == true
    end

    @testset "set_flow_partitions_num!" begin
        data = WaterModels.parse_file("../examples/data/epanet/van_zyl.inp")
        WaterModels.set_flow_partitions_num!(data, 5)
        @test length(data["pipe"]["1"]["flow_partition"]) == 5
        @test length(data["pump"]["1"]["flow_partition"]) == 5
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

    @testset "set_warm_start! (multinetwork)" begin
        network = WaterModels.parse_file("../examples/data/epanet/van_zyl.inp")
        network_mn = WaterModels.make_multinetwork(network)
        WaterModels.set_warm_start!(network_mn)
        @test haskey(network_mn["nw"]["1"]["pipe"]["1"], "q_start")
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
        base_flow = WaterModels._calc_flow_per_unit_transform(network_data)(0.02777)
        @test isapprox(network_data["time_series"]["demand"]["2"]["flow_nominal"][1], base_flow, rtol=1.0e-4)
        @test isapprox(network_data["time_series"]["demand"]["2"]["flow_nominal"][2], 0.5*base_flow, rtol=1.0e-4)
        @test isapprox(network_data["time_series"]["demand"]["2"]["flow_nominal"][3], 0.25*base_flow, rtol=1.0e-4)

        ts_length = network_data["time_series"]["num_steps"]
        mn_data = WaterModels.make_multinetwork(network_data)

        @test _IM.ismultinetwork(mn_data)
        @test !haskey(mn_data, "time_series")
        @test length(mn_data["nw"]) == ts_length
        @test isapprox(mn_data["nw"]["1"]["demand"]["2"]["flow_nominal"], base_flow, rtol=1.0e-4)
        @test isapprox(mn_data["nw"]["2"]["demand"]["2"]["flow_nominal"], 0.5*base_flow, rtol=1.0e-4)
        @test isapprox(mn_data["nw"]["3"]["demand"]["2"]["flow_nominal"], 0.25*base_flow, rtol=1.0e-4)
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
        WaterModels.correct_network_data!(network_data)
        WaterModels.convert_short_pipes!(network_data)
        @test haskey(network_data["short_pipe"], "1")

        network_data = WaterModels.parse_epanet("../test/data/epanet/snapshot/short-pipe-valve-lps.inp")
        WaterModels.epanet_to_watermodels!(network_data)
        WaterModels.correct_network_data!(network_data)
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

    @testset "propagate_topology_status!" begin
        @testset "inactive demand by propagation" begin
            network = WaterModels.parse_file("../test/data/epanet/snapshot/pipe-hw-lps.inp")
            network["node"]["2"]["status"] = WaterModels.STATUS_INACTIVE
            WaterModels.propagate_topology_status!(network)
            @test network["demand"]["2"]["status"] === WaterModels.STATUS_INACTIVE
        end

        @testset "inactive reservoir by propagation" begin
            network = WaterModels.parse_file("../test/data/epanet/snapshot/pipe-hw-lps.inp")
            network["node"]["1"]["status"] = WaterModels.STATUS_INACTIVE
            WaterModels.propagate_topology_status!(network)
            @test network["reservoir"]["1"]["status"] === WaterModels.STATUS_INACTIVE
        end

        @testset "inactive tank by propagation" begin
            network = WaterModels.parse_file("../test/data/epanet/snapshot/tank-hw-lps.inp")
            network["node"]["1"]["status"] = WaterModels.STATUS_INACTIVE
            WaterModels.propagate_topology_status!(network)
            @test network["tank"]["1"]["status"] === WaterModels.STATUS_INACTIVE
        end

        @testset "inactive pipe by propagation" begin
            network = WaterModels.parse_file("../test/data/epanet/snapshot/pipe-hw-lps.inp")
            network["node"]["1"]["status"] = WaterModels.STATUS_INACTIVE
            WaterModels.propagate_topology_status!(network)
            @test network["pipe"]["1"]["status"] === WaterModels.STATUS_INACTIVE
        end

        @testset "inactive design pipe by propagation" begin
            network = WaterModels.parse_file("../test/data/json/shamir.json")
            network["node"]["2"]["status"] = WaterModels.STATUS_INACTIVE
            
            @test WaterModels.propagate_topology_status!(network) == true
            @test network["des_pipe"]["1"]["status"] === WaterModels.STATUS_INACTIVE
        end

        @testset "inactive valve by propagation" begin
            network = WaterModels.parse_file("../test/data/epanet/snapshot/shutoff_valve-hw-lps.inp")
            network["node"]["3"]["status"] = WaterModels.STATUS_INACTIVE
            
            @test WaterModels.propagate_topology_status!(network) == true
            @test network["valve"]["1"]["status"] === WaterModels.STATUS_INACTIVE
        end

        @testset "inactive regulator by propagation" begin
            network = WaterModels.parse_file("../test/data/epanet/snapshot/prv-hw-lps.inp")
            network["node"]["3"]["status"] = WaterModels.STATUS_INACTIVE
            
            @test WaterModels.propagate_topology_status!(network) == true
            @test network["regulator"]["2"]["status"] === WaterModels.STATUS_INACTIVE
        end

        @testset "inactive short pipe by propagation" begin
            network = WaterModels.parse_file("../test/data/epanet/snapshot/short-pipe-lps.inp")
            WaterModels.convert_short_pipes!(network)
            network["node"]["2"]["status"] = WaterModels.STATUS_INACTIVE
            
            @test WaterModels.propagate_topology_status!(network) == true
            @test network["short_pipe"]["1"]["status"] === WaterModels.STATUS_INACTIVE
        end

        @testset "inactive pump by propagation" begin
            network = WaterModels.parse_file("../test/data/epanet/snapshot/pump-hw-lps.inp")
            network["node"]["2"]["status"] = WaterModels.STATUS_INACTIVE

            @test WaterModels.propagate_topology_status!(network) == true
            @test network["pump"]["1"]["status"] === WaterModels.STATUS_INACTIVE
        end

        @testset "multinetwork variant" begin
            network = WaterModels.parse_file("../test/data/epanet/multinetwork/pipe-hw-lps.inp")
            network_mn = WaterModels.make_multinetwork(network)
            network_mn["nw"]["1"]["node"]["1"]["status"] = WaterModels.STATUS_INACTIVE
            
            @test WaterModels.propagate_topology_status!(network_mn) == true
            @test network_mn["nw"]["1"]["pipe"]["1"]["status"] === WaterModels.STATUS_INACTIVE
            @test network_mn["nw"]["2"]["pipe"]["1"]["status"] === WaterModels.STATUS_ACTIVE
        end
    end

    @testset "simplify_network!" begin
        network = WaterModels.parse_file("../test/data/epanet/multinetwork/pipe-hw-lps.inp")
        network_mn = WaterModels.make_multinetwork(network)
        network_mn["nw"]["1"]["node"]["1"]["status"] = WaterModels.STATUS_INACTIVE
        WaterModels.simplify_network!(network_mn) == true
        
        @test network_mn["nw"]["1"]["pipe"]["1"]["status"] === WaterModels.STATUS_INACTIVE
        @test network_mn["nw"]["2"]["pipe"]["1"]["status"] === WaterModels.STATUS_ACTIVE
    end
end
