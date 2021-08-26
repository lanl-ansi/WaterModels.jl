@testset "src/core/node.jl" begin
    @testset "_relax_nodes! (single network with `time_series`)" begin
        data = WaterModels.parse_file("../test/data/epanet/multinetwork/reservoir-hw-lps.inp")
        WaterModels._relax_nodes!(data)
        @test isapprox(data["node"]["1"]["head_min"] * data["base_head"], 10.0)
        @test isapprox(data["node"]["1"]["head_max"] * data["base_head"], 20.0)
    end
end
