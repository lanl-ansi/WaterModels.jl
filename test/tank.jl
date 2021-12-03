@testset "src/core/tank.jl" begin
    @testset "make_all_nondispatchable! (multinetwork)" begin
        data = WaterModels.parse_file("../test/data/epanet/multinetwork/tank-hw-lps.inp")
        mn_data = WaterModels.make_multinetwork(data)

        mn_data["nw"]["1"]["tank"]["1"]["dispatchable"] = true
        mn_data["nw"]["3"]["tank"]["1"]["dispatchable"] = true

        WaterModels.make_all_nondispatchable!(mn_data)
        @test mn_data["nw"]["1"]["tank"]["1"]["dispatchable"] == false
        @test mn_data["nw"]["3"]["tank"]["1"]["dispatchable"] == false
    end

    @testset "set_tank_bounds_from_time_series! (single network)" begin
        data = WaterModels.parse_file("../test/data/epanet/multinetwork/tank-hw-lps.inp")
        WaterModels.set_tank_bounds_from_time_series!(data)
        @test data["tank"]["1"]["dispatchable"] == true 
        @test data["tank"]["1"]["dispatchable"] == true
    end

    @testset "set_tank_bounds_from_time_series! (multinetwork)" begin
        data = WaterModels.parse_file("../test/data/epanet/multinetwork/tank-hw-lps.inp")
        mn_data = WaterModels.make_multinetwork(data)
        @test_throws AssertionError WaterModels.set_tank_bounds_from_time_series!(mn_data)
    end

    @testset "make_tank_start_dispatchable! (single network)" begin
        data = WaterModels.parse_file("../test/data/epanet/snapshot/tank-hw-lps.inp")
        WaterModels.make_tank_start_dispatchable!(data)
        @test data["tank"]["1"]["dispatchable"] == true
    end

    @testset "make_tank_start_dispatchable! (multinetwork)" begin
        data = WaterModels.parse_file("../test/data/epanet/multinetwork/tank-hw-lps.inp")
        mn_data = WaterModels.make_multinetwork(data)
        WaterModels.make_tank_start_dispatchable!(mn_data)
        @test mn_data["nw"]["1"]["tank"]["1"]["dispatchable"] == true 
        @test mn_data["nw"]["3"]["tank"]["1"]["dispatchable"] == false
    end
end
