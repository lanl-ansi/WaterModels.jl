@testset "src/io/common.jl" begin
    @testset "parse_file (.inp)" begin
        balerma_data = WaterModels.parse_file("../test/data/epanet/balerma.inp")
        balerma_name = "Balerma Network"
        @test balerma_data["name"] == balerma_name

        klmod_data = WaterModels.parse_file("../test/data/epanet/klmod.inp")
        klmod_name = "Global Water Full network - Peak Day (Avg * 1.9)"
        @test klmod_data["name"] == klmod_name

        richmond_skeleton_data = WaterModels.parse_file("../test/data/epanet/richmond-skeleton.inp")
        richmond_skeleton_name = "Richmond Skeleton Water Supply System"
        @test richmond_skeleton_data["name"] == richmond_skeleton_name

        richmond_skeleton_sp_data = WaterModels.parse_file("../test/data/epanet/richmond-skeleton-sp.inp")
        richmond_skeleton_sp_name = "Richmond Skeleton Water Supply System (First Time Point)"
        @test richmond_skeleton_sp_data["name"] == richmond_skeleton_sp_name

        shamir_data = WaterModels.parse_file("../test/data/epanet/shamir.inp")
        shamir_name = "shamir -- Bragalli, D'Ambrosio, Lee, Lodi, Toth (2008)"
        @test shamir_data["name"] == shamir_name

        shamir_ts_data = WaterModels.parse_file("../test/data/epanet/shamir-ts.inp")
        shamir_ts_name = "shamir (time series) -- Bragalli, D'Ambrosio, Lee, Lodi, Toth (2008)"
        @test shamir_ts_data["name"] == shamir_ts_name
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
