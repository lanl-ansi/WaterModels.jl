@testset "src/core/demand.jl" begin
    @testset "make_all_nondispatchable! (multinetwork)" begin
        data = WaterModels.parse_file("../test/data/epanet/multinetwork/pipe-hw-lps.inp")
        mn_data = WaterModels.make_multinetwork(data)
        WaterModels.make_all_nondispatchable!(mn_data)
        @test mn_data["nw"]["1"]["demand"]["2"]["dispatchable"] == false
        @test mn_data["nw"]["3"]["demand"]["2"]["dispatchable"] == false
    end

    @testset "set_demand_bounds_from_time_series! (single network with `time_series`)" begin
        data = WaterModels.parse_file("../test/data/epanet/multinetwork/pipe-hw-lps.inp")
        WaterModels.set_demand_bounds_from_time_series!(data)
        @test data["demand"]["2"]["dispatchable"] == true
    end

    @testset "set_demand_bounds_from_time_series! (multinetwork)" begin
        data = WaterModels.parse_file("../test/data/epanet/multinetwork/pipe-hw-lps.inp")
        mn_data = WaterModels.make_multinetwork(data)
        @test_throws AssertionError WaterModels.set_demand_bounds_from_time_series!(mn_data)
    end
end
