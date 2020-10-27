@testset "src/core/node.jl" begin
    @testset "_relax_nodes! (single network with `time_series`)" begin
        data = WaterModels.parse_file("../test/data/epanet/multinetwork/reservoir-hw-lps.inp")
        WaterModels._relax_nodes!(data)
        @test data["node"]["1"]["h_min"] == 10.0
        @test data["node"]["1"]["h_max"] == 20.0
    end
end
