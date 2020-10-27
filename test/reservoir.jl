@testset "src/core/reservoir.jl" begin
    @testset "_relax_reservoirs! (single network)" begin
        data = WaterModels.parse_file("../test/data/epanet/snapshot/pipe-hw-lps.inp")
        WaterModels._relax_reservoirs!(data)
        @test data["reservoir"]["1"]["dispatchable"] == true
    end

    @testset "_relax_reservoirs! (multinetwork)" begin
        data = WaterModels.parse_file("../test/data/epanet/multinetwork/pipe-hw-lps.inp")
        mn_data = WaterModels.make_multinetwork(data)
        WaterModels._relax_reservoirs!(mn_data)
        @test mn_data["nw"]["1"]["reservoir"]["1"]["dispatchable"] == true 
        @test mn_data["nw"]["3"]["reservoir"]["1"]["dispatchable"] == true
    end
end
