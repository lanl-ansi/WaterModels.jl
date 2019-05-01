@testset "src/io/common.jl" begin
    @testset "parse_file (.inp)" begin
        shamir_data = WaterModels.parse_file("../test/data/epanet/shamir.inp")
        shamir_title = "shamir -- Bragalli, D'Ambrosio, Lee, Lodi, Toth (2008)"
        @test shamir_data["title"] == lowercase(shamir_title)

        balerma_data = WaterModels.parse_file("../test/data/epanet/balerma.inp")
        balerma_title = "Balerma Network"
        @test balerma_data["title"] == lowercase(balerma_title)
    end

    @testset "parse_file (.json)" begin
        data = WaterModels.parse_file("../test/data/json/shamir.json")
        pipe_1_max_velocity = data["pipes"]["1"]["maximumVelocity"]
        @test pipe_1_max_velocity == 2.0
    end

    @testset "parse_file (invalid extension)" begin
        path = "../test/data/json/shamir.data"
        @test_throws ErrorException WaterModels.parse_file(path)
    end

    @testset "parse_json" begin
        data = WaterModels.parse_json("../test/data/json/shamir.json")
        pipe_1_max_velocity = data["pipes"]["1"]["maximumVelocity"]
        @test pipe_1_max_velocity == 2.0
    end
end
