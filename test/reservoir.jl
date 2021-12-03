@testset "src/core/reservoir.jl" begin
    @testset "make_all_nondispatchable! (multinetwork)" begin
        data = WaterModels.parse_file("../test/data/epanet/multinetwork/reservoir-hw-lps.inp")
        mn_data = WaterModels.make_multinetwork(data)

        mn_data["nw"]["1"]["reservoir"]["1"]["dispatchable"] = true
        mn_data["nw"]["3"]["reservoir"]["1"]["dispatchable"] = true

        WaterModels.make_all_nondispatchable!(mn_data)
        @test mn_data["nw"]["1"]["reservoir"]["1"]["dispatchable"] == false
        @test mn_data["nw"]["3"]["reservoir"]["1"]["dispatchable"] == false
    end

    @testset "set_reservoir_bounds_from_time_series! (single network)" begin
        data = WaterModels.parse_file("../test/data/epanet/multinetwork/reservoir-hw-lps.inp")
        WaterModels.set_reservoir_bounds_from_time_series!(data)
        @test data["reservoir"]["1"]["dispatchable"] == true
    end

    @testset "set_reservoir_bounds_from_time_series! (multinetwork)" begin
        data = WaterModels.parse_file("../test/data/epanet/multinetwork/reservoir-hw-lps.inp")
        mn_data = WaterModels.make_multinetwork(data)
        @test_throws AssertionError WaterModels.set_reservoir_bounds_from_time_series!(mn_data)
    end
end
