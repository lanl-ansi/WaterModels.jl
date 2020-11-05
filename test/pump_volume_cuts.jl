@testset "src/util/pump_volume_cuts.jl" begin
    @testset "_add_pump_volume_cuts!" begin
        data = WaterModels.parse_file("../test/data/epanet/multinetwork/owf-hw-lps.inp")
        mn_data = WaterModels.make_multinetwork(data)
        wm = instantiate_model(mn_data, LRDWaterModel, build_mn_wf)
        WaterModels._add_pump_volume_cuts!(wm) # Add the cuts to the model.
        result = WaterModels.optimize_model!(wm, optimizer=cbc)
        @test result["termination_status"] == OPTIMAL
    end
end
