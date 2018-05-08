@testset "test minlp feasibility problems" begin
    @testset "balerma" begin
        result = build_generic_model("../test/data/epanet/balerma.inp", GenericWaterModel)
    end
end
