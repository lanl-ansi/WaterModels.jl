@testset "src/core/function.jl" begin
    @testset "_get_alpha_min_1 (with non-matching exponents)" begin
        data = WaterModels.parse_file("../test/data/epanet/multinetwork/pipe-hw-lps.inp")
        mn_data = WaterModels.make_multinetwork(data)
        wm = instantiate_model(mn_data, NCWaterModel, build_mn_wf)
        ref(wm, 1)[:alpha] = NaN # Change one of the exponents.
        @test_throws ErrorException WaterModels._get_alpha_min_1(wm)
    end
end
