@testset "test minlp expansion problems" begin
    @testset "hanoi" begin
        solution = run_expansion("../test/data/epanet/hanoi.inp", GenericWaterModel{StandardMINLPForm}, solver)
        @test solution["status"] == :LocalOptimal
    end
end
