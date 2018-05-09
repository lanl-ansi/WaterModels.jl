@testset "test minlp feasibility problems" begin
    @testset "balerma" begin
        model = build_generic_model("../test/data/epanet/balerma.inp", GenericWaterModel{StandardMINLPForm})
        @test true
    end
end
