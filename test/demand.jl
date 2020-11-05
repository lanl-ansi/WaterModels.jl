@testset "src/core/demand.jl" begin
    @testset "_relax_demand! (has keys)" begin
        data = WaterModels.parse_file("../test/data/epanet/snapshot/pipe-hw-lps.inp")
        data["demand"]["2"]["demand_min"], data["demand"]["2"]["demand_max"] = 0.0, 1.0
        WaterModels._relax_demand!(data["demand"]["2"])
        @test data["demand"]["2"]["dispatchable"] == true
    end

    @testset "_relax_demand! (does not have keys)" begin
        data = WaterModels.parse_file("../test/data/epanet/snapshot/pipe-hw-lps.inp")
        WaterModels._relax_demand!(data["demand"]["2"])
        @test data["demand"]["2"]["dispatchable"] == false
    end

    @testset "_relax_demands! (single network)" begin
        data = WaterModels.parse_file("../test/data/epanet/snapshot/pipe-hw-lps.inp")
        WaterModels._relax_demands!(data)
        @test data["demand"]["2"]["dispatchable"] == false
    end

    @testset "_relax_demands! (single network with `time_series`)" begin
        data = WaterModels.parse_file("../test/data/epanet/multinetwork/pipe-hw-lps.inp")
        WaterModels._relax_demands!(data)
        @test data["demand"]["2"]["dispatchable"] == true
    end

    @testset "_relax_demands! (multinetwork)" begin
        data = WaterModels.parse_file("../test/data/epanet/multinetwork/pipe-hw-lps.inp")
        mn_data = WaterModels.make_multinetwork(data)
        WaterModels._relax_demands!(mn_data)
        @test mn_data["nw"]["1"]["demand"]["2"]["dispatchable"] == false
        @test mn_data["nw"]["3"]["demand"]["2"]["dispatchable"] == false
    end
end
