@testset "src/io/common.jl" begin
    @testset "parse_file (.inp)" begin
        balerma_data = WaterModels.parse_file("../test/data/epanet/balerma.inp")
        balerma_title = "Balerma Network"
        @test balerma_data["title"] == balerma_title

        klmod_data = WaterModels.parse_file("../test/data/epanet/klmod.inp")
        klmod_title = "Global Water Full network - Peak Day (Avg * 1.9)"
        @test klmod_data["title"] == klmod_title

        richmond_skeleton_data = WaterModels.parse_file("../test/data/epanet/richmond-skeleton.inp")
        richmond_skeleton_title = "Richmond Skeleton Water Supply System"
        @test richmond_skeleton_data["title"] == richmond_skeleton_title

        shamir_data = WaterModels.parse_file("../test/data/epanet/shamir.inp")
        shamir_title = "shamir -- Bragalli, D'Ambrosio, Lee, Lodi, Toth (2008)"
        @test shamir_data["title"] == shamir_title
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
