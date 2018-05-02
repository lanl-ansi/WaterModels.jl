@testset "test data summary" begin
    @testset "balerma from EPANET file" begin
        data = WaterModels.parse_file("../test/data/epanet/balerma.inp")
    end

    @testset "foss_poly_1 from EPANET file" begin
        data = WaterModels.parse_file("../test/data/epanet/foss_poly_1.inp")
    end

    @testset "hanoi_extended from EPANET file" begin
        data = WaterModels.parse_file("../test/data/epanet/hanoi_extended.inp")
    end

    @testset "hanoi from EPANET file" begin
        data = WaterModels.parse_file("../test/data/epanet/hanoi.inp")
    end

    @testset "klmod from EPANET file" begin
        data = WaterModels.parse_file("../test/data/epanet/klmod.inp")
    end

    @testset "rural from EPANET file" begin
        data = WaterModels.parse_file("../test/data/epanet/rural.inp")
    end

    @testset "zj from EPANET file" begin
        data = WaterModels.parse_file("../test/data/epanet/zj.inp")
    end
end
