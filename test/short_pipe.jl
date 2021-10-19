@testset "src/core/short_pipe.jl" begin
    @testset "correct_short_pipes!" begin
        network = WaterModels.parse_file("../test/data/epanet/snapshot/short-pipe-lps.inp")
        WaterModels.convert_short_pipes!(network)
        WaterModels.correct_short_pipes!(network)
    end

    @testset "set_short_pipe_warm_start!" begin
        network = WaterModels.parse_file("../test/data/epanet/snapshot/short-pipe-lps.inp")
        WaterModels.convert_short_pipes!(network)
        WaterModels.correct_short_pipes!(network)
        WaterModels.set_short_pipe_warm_start!(network)
    end
end