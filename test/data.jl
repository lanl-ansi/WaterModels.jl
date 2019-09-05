@testset "src/core/data.jl" begin
    network_path_shamir = "../test/data/epanet/shamir-ts.inp"
    network_path_richmond = "../test/data/epanet/richmond-skeleton.inp"

    @testset "make_multinetwork shamir" begin
        # This should silence everything except error messages.
        network_data = WaterModels.parse_file(network_path_shamir)
        @test !InfrastructureModels.ismultinetwork(network_data)
        @test haskey(network_data, "time_series")

        ts_length = network_data["time_series"]["num_steps"]

        mn_data = WaterModels.make_multinetwork(network_data)
        @test InfrastructureModels.ismultinetwork(mn_data)
        @test !haskey(mn_data, "time_series")
        @test length(mn_data["nw"]) == ts_length
    end

    @testset "make_multinetwork richmond" begin
        # This should silence everything except error messages.
        network_data = WaterModels.parse_file(network_path_richmond)
        @test !InfrastructureModels.ismultinetwork(network_data)
        @test haskey(network_data, "time_series")

        ts_length = network_data["time_series"]["num_steps"]

        mn_data = WaterModels.make_multinetwork(network_data)
        @test InfrastructureModels.ismultinetwork(mn_data)
        @test !haskey(mn_data, "time_series")
        @test length(mn_data["nw"]) == ts_length

        junction_ids = keys(network_data["time_series"]["junction"])
        for step_index in 1:ts_length
            ts_demand = sum(network_data["time_series"]["junction"][jid]["demand"][step_index] for jid in junction_ids)
            mn_demand = sum(mn_data["nw"]["$(step_index)"]["junction"][jid]["demand"] for jid in junction_ids)
            @test isapprox(ts_demand, mn_demand)
        end
    end

    @testset "load_timepoint! richmond" begin
        # This should silence everything except error messages.
        network_data = WaterModels.parse_file(network_path_richmond)
        @test !InfrastructureModels.ismultinetwork(network_data)
        @test haskey(network_data, "time_series")

        ts_length = network_data["time_series"]["num_steps"]

        junction_ids = keys(network_data["time_series"]["junction"])
        for step_index in 1:ts_length
            InfrastructureModels.load_timepoint!(network_data, step_index)
            ts_demand = sum(network_data["time_series"]["junction"][jid]["demand"][step_index] for jid in junction_ids)
            demand = sum(network_data["junction"][jid]["demand"] for jid in junction_ids)
            @test isapprox(ts_demand, demand)
        end
    end

end
