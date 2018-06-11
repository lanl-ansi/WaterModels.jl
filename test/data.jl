@testset "test data summary" begin
    @testset "balerma from EPANET file" begin
        data = WaterModels.parse_file("../test/data/epanet/balerma.inp")
        @test data["title"] == lowercase("Balerma Network")
    end

    @testset "d-town from EPANET file" begin
        data = WaterModels.parse_file("../test/data/epanet/d-town.inp")
        @test data["title"] == lowercase("")
    end

    @testset "foss_poly_1 from EPANET file" begin
        data = WaterModels.parse_file("../test/data/epanet/foss_poly_1.inp")
        @test data["title"] == lowercase("foss_poly_1 -- Bragalli, D'Ambrosio, Lee, Lodi, Toth (2008)")
    end

    @testset "hanoi_extended from EPANET file" begin
        data = WaterModels.parse_file("../test/data/epanet/hanoi_extended.inp")
        @test data["title"] == lowercase("")
    end

    @testset "hanoi from EPANET file" begin
        data = WaterModels.parse_file("../test/data/epanet/hanoi.inp")
        @test data["title"] == lowercase("")
    end

    @testset "klmod from EPANET file" begin
        data = WaterModels.parse_file("../test/data/epanet/klmod.inp")
        @test data["title"] == lowercase("Global Water Full network - Peak Day (Avg * 1.9)")
    end

    @testset "rural from EPANET file" begin
        data = WaterModels.parse_file("../test/data/epanet/rural.inp")
        @test data["title"] == lowercase("")
    end

    @testset "tasseff from EPANET file" begin
        data = WaterModels.parse_file("../test/data/epanet/tasseff.inp")
        @test data["title"] == lowercase("")
    end

    @testset "zj from EPANET file" begin
        data = WaterModels.parse_file("../test/data/epanet/zj.inp")
        @test data["title"] == lowercase("")
    end
end
