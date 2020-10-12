@testset "src/io/common.jl" begin
    @testset "parse_file (.inp)" begin
        lps = WaterModels.parse_file("../test/data/epanet/snapshot/pipe-hw-lps.inp")
        gpm = WaterModels.parse_file("../test/data/epanet/snapshot/pipe-hw-gpm.inp")

        @test isapprox(lps["demand"]["2"]["flow_rate"], gpm["demand"]["2"]["flow_rate"], rtol=1.0e-4)
        @test isapprox(lps["node"]["1"]["elevation"], gpm["node"]["1"]["elevation"], rtol=1.0e-4)

        balerma_data = WaterModels.parse_file("../test/data/epanet/balerma.inp")
        balerma_name = "Balerma Network"
        @test balerma_data["name"] == balerma_name

        klmod_data = WaterModels.parse_file("../test/data/epanet/klmod.inp")
        klmod_name = "Global Water Full network - Peak Day (Avg * 1.9)"
        @test klmod_data["name"] == klmod_name

        shamir_data = WaterModels.parse_file("../examples/data/epanet/shamir.inp")
        shamir_name = "Shamir (Two-loop) Water Network"
        @test shamir_data["name"] == shamir_name

        shamir_ts_data = WaterModels.parse_file("../test/data/epanet/shamir-ts.inp")
        shamir_ts_name = "shamir (time series) -- Bragalli, D'Ambrosio, Lee, Lodi, Toth (2008)"
        @test shamir_ts_data["name"] == shamir_ts_name
    end

    @testset "parse_file (.json)" begin
        data = WaterModels.parse_file("../examples/data/json/shamir.json")
        pipe_1_max_velocity = data["pipe"]["1"]["maximumVelocity"]
        @test pipe_1_max_velocity == 2.0
    end

    @testset "parse_file (invalid extension)" begin
        path = "../examples/data/json/shamir.data"
        @test_throws ErrorException WaterModels.parse_file(path)
    end

    @testset "parse_json" begin
        data = WaterModels.parse_json("../examples/data/json/shamir.json")
        pipe_1_max_velocity = data["pipe"]["1"]["maximumVelocity"]
        @test pipe_1_max_velocity == 2.0
    end
end
