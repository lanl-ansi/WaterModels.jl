# Iterate over all possible WaterModels formulations.
for formulation in [NCWaterModel, NCDWaterModel, CRDWaterModel, LAWaterModel, LRDWaterModel, PWLRDWaterModel]
    # Set a generic extensions dictionary for testing purposes.
    ext = Dict(:pipe_breakpoints => 3, :pump_breakpoints => 3)

    @testset "Network Design Problems: $(formulation)" begin
        network = parse_file("../test/data/json/shamir.json")

        @testset "Shamir Network Design (Reduced): $(formulation)" begin
            wm = instantiate_model(network, formulation, build_des; ext = ext)
            result = WaterModels.optimize_model!(wm, optimizer=_choose_solver(wm, ipopt, cbc))
            @test _is_valid_status(result["termination_status"])
            @test isapprox(result["objective"], 1.36e6, rtol=1.0e-4)
        end
    end
end

@testset "solve_des" begin
    network = WaterModels.parse_file("../test/data/json/shamir.json")
    result = WaterModels.solve_des(network, LRDWaterModel, cbc)
    result = WaterModels.run_des(network, LRDWaterModel, cbc)
    @test result["termination_status"] == OPTIMAL
end

@testset "solve_des (per-unit transformation)" begin
    # Note that the per-unit transformation is currently experimental.
    network = WaterModels.parse_file("../test/data/json/shamir.json")
    make_per_unit!(network) # This function is currently experimental.
    result = WaterModels.solve_des(network, LRDWaterModel, cbc)
    result = WaterModels.run_des(network, LRDWaterModel, cbc)
    @test result["termination_status"] == OPTIMAL
end