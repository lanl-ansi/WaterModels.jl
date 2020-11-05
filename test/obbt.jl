@testset "src/util/obbt.jl" begin
    @testset "run_obbt_owf! (pipe)" begin
        data = WaterModels.parse_file("../test/data/epanet/multinetwork/pipe-hw-lps.inp")
        run_obbt_owf!(data, cbc; model_type = LRDWaterModel, max_iter = 1)
        @test haskey(data["pipe"]["1"], "q_min")
        @test haskey(data["pipe"]["1"], "q_max")
    end

    @testset "run_obbt_owf! (pump)" begin
        data = WaterModels.parse_file("../test/data/epanet/multinetwork/pump-hw-lps.inp")
        run_obbt_owf!(data, cbc; model_type = LRDWaterModel, max_iter = 1)
        @test haskey(data["pump"]["1"], "q_max")
    end

    @testset "run_obbt_owf! (regulator)" begin
        data = WaterModels.parse_file("../test/data/epanet/multinetwork/prv-hw-lps.inp")
        run_obbt_owf!(data, cbc; model_type = LRDWaterModel, max_iter = 1)
        @test haskey(data["regulator"]["2"], "q_max")
    end

    @testset "run_obbt_owf! (short pipe)" begin
        data = WaterModels.parse_file("../test/data/epanet/multinetwork/short-pipe-lps.inp")
        run_obbt_owf!(data, cbc; model_type = LRDWaterModel, max_iter = 1)
        @test haskey(data["short_pipe"]["1"], "q_min")
        @test haskey(data["short_pipe"]["1"], "q_max")
    end

    @testset "run_obbt_owf! (valve)" begin
        data = WaterModels.parse_file("../test/data/epanet/multinetwork/shutoff_valve-hw-lps.inp")
        run_obbt_owf!(data, cbc; model_type = LRDWaterModel, max_iter = 1)
        @test haskey(data["valve"]["1"], "q_min")
        @test haskey(data["valve"]["1"], "q_max")
    end
end
