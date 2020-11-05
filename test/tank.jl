@testset "src/core/tank.jl" begin
    @testset "_relax_tanks! (single network)" begin
        data = WaterModels.parse_file("../test/data/epanet/snapshot/tank-hw-lps.inp")
        WaterModels._relax_tanks!(data)
        @test data["tank"]["1"]["dispatchable"] == true
    end

    @testset "_relax_tanks! (multinetwork)" begin
        data = WaterModels.parse_file("../test/data/epanet/multinetwork/tank-hw-lps.inp")
        mn_data = WaterModels.make_multinetwork(data)
        WaterModels._relax_tanks!(mn_data)
        @test mn_data["nw"]["1"]["tank"]["1"]["dispatchable"] == true 
        @test mn_data["nw"]["3"]["tank"]["1"]["dispatchable"] == true
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
