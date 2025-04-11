# Iterate over all possible WaterModels formulations.
for formulation in [NCWaterModel, NCDWaterModel, CRDWaterModel, LAWaterModel, LRDWaterModel, PWLRDWaterModel]
    @testset "Network Design Problems: $(formulation)" begin
        network = parse_file("../test/data/json/shamir.json")
        set_flow_partitions_si!(network, 10.0, 1.0e-4)

        @testset "Shamir Network Design (Reduced): $(formulation)" begin
            wm = instantiate_model(network, formulation, build_des)
            result = WaterModels.optimize_model!(wm, optimizer = _choose_solver(wm, nlp_solver, milp_solver))
            if _is_valid_status(result["termination_status"])
                @test isapprox(result["objective"], 1.36e6, rtol=1.0e-4)
            else
                # FIXME(odow): this test is broken ubuntu and windows in Julia
                # v1.11 but not v1.6, and not on macOS.
                @test_broken isapprox(result["objective"], 1.36e6, rtol=1.0e-4)
            end
        end
    end
end

@testset "solve_des" begin
    network = WaterModels.parse_file("../test/data/json/shamir.json")
    set_flow_partitions_si!(network, 10.0, 1.0e-4)

    result = WaterModels.solve_des(network, LRDWaterModel, milp_solver)
    result = WaterModels.run_des(network, LRDWaterModel, milp_solver)
    @test result["termination_status"] == OPTIMAL
end
