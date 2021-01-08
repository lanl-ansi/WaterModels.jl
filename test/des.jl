# Iterate over all possible WaterModels formulations.
for formulation in [NCWaterModel, NCDWaterModel, CRDWaterModel, LAWaterModel, LRDWaterModel, PWLRDWaterModel]
    # Set a generic extensions dictionary for testing purposes.
    ext = Dict(:pipe_breakpoints => 3, :pump_breakpoints => 3)

    @testset "Network Design Problems: $(formulation)" begin
        network = parse_json("../test/data/json/shamir.json")

        @testset "Shamir Network Design (Reduced): $(formulation)" begin
            wm = instantiate_model(network, formulation, build_des; ext = ext)
            result = WaterModels.optimize_model!(wm, optimizer=_make_juniper(wm, ipopt))
            @test result["termination_status"] == LOCALLY_SOLVED
            @test isapprox(result["objective"], 1.36e6, rtol=1.0e-4)
        end
    end
end