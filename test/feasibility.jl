@testset "test minlp feasibility problems" begin
    @testset "balerma" begin
        model = build_generic_model("../test/data/epanet/balerma.inp", GenericWaterModel{StandardMINLPForm})
        @test true
    end

    @testset "foss_poly_1" begin
        model = build_generic_model("../test/data/epanet/foss_poly_1.inp", GenericWaterModel{StandardMINLPForm})
        @test true
    end

    @testset "hanoi_extended" begin
        model = build_generic_model("../test/data/epanet/hanoi_extended.inp", GenericWaterModel{StandardMINLPForm})
        @test true
    end

    @testset "hanoi" begin
        model = build_generic_model("../test/data/epanet/hanoi.inp", GenericWaterModel{StandardMINLPForm})
        @test true
    end

    @testset "klmod" begin
        model = build_generic_model("../test/data/epanet/klmod.inp", GenericWaterModel{StandardMINLPForm})
        @test true
    end

    @testset "rural" begin
        model = build_generic_model("../test/data/epanet/rural.inp", GenericWaterModel{StandardMINLPForm})
        @test true
    end

    @testset "zj" begin
        model = build_generic_model("../test/data/epanet/zj.inp", GenericWaterModel{StandardMINLPForm})
        @test true
    end
end
