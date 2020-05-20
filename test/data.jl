@testset "src/core/data.jl" begin
    @testset "make_multinetwork shamir" begin
        network_data = WaterModels.parse_file("../test/data/epanet/shamir-ts.inp")

        @test !InfrastructureModels.ismultinetwork(network_data)
        @test haskey(network_data, "time_series")
        @test isapprox(network_data["time_series"]["junction"]["2"]["demand"][1], 0.02777, rtol=1.0e-4)
        @test isapprox(network_data["time_series"]["junction"]["2"]["demand"][2], 0.5*0.02777, rtol=1.0e-4)
        @test isapprox(network_data["time_series"]["junction"]["2"]["demand"][3], 0.25*0.02777, rtol=1.0e-4)

        ts_length = network_data["time_series"]["num_steps"]
        mn_data = WaterModels.make_multinetwork(network_data)

        @test InfrastructureModels.ismultinetwork(mn_data)
        @test !haskey(mn_data, "time_series")
        @test length(mn_data["nw"]) == ts_length
        @test isapprox(mn_data["nw"]["1"]["junction"]["2"]["demand"], 0.02777, rtol=1.0e-4)
        @test isapprox(mn_data["nw"]["2"]["junction"]["2"]["demand"], 0.5*0.02777, rtol=1.0e-4)
        @test isapprox(mn_data["nw"]["3"]["junction"]["2"]["demand"], 0.25*0.02777, rtol=1.0e-4)
    end
end
