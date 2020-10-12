@testset "src/core/data.jl" begin
    @testset "make_multinetwork shamir" begin
        network_data = WaterModels.parse_file("../test/data/epanet/shamir-ts.inp")

        @test !_IM.ismultinetwork(network_data)
        @test haskey(network_data, "time_series")
        @test isapprox(network_data["time_series"]["demand"]["2"]["flow_rate"][1], 0.02777, rtol=1.0e-4)
        @test isapprox(network_data["time_series"]["demand"]["2"]["flow_rate"][2], 0.5*0.02777, rtol=1.0e-4)
        @test isapprox(network_data["time_series"]["demand"]["2"]["flow_rate"][3], 0.25*0.02777, rtol=1.0e-4)

        ts_length = network_data["time_series"]["num_steps"]
        mn_data = WaterModels.make_multinetwork(network_data)

        @test _IM.ismultinetwork(mn_data)
        @test !haskey(mn_data, "time_series")
        @test length(mn_data["nw"]) == ts_length
        @test isapprox(mn_data["nw"]["1"]["demand"]["2"]["flow_rate"], 0.02777, rtol=1.0e-4)
        @test isapprox(mn_data["nw"]["2"]["demand"]["2"]["flow_rate"], 0.5*0.02777, rtol=1.0e-4)
        @test isapprox(mn_data["nw"]["3"]["demand"]["2"]["flow_rate"], 0.25*0.02777, rtol=1.0e-4)
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
end
