@testset "src/core/pump.jl" begin
    for head_curve_form in [PUMP_QUADRATIC, PUMP_BEST_EFFICIENCY_POINT, PUMP_EPANET, PUMP_LINEAR_POWER]
        @testset "pump with efficiency curve: $(head_curve_form)" begin
            # Read in the initial network data and ensure an efficiency curve exists.
            network_path = "../test/data/epanet/snapshot/pump-hw-lps-eff-curve.inp"
            network = WaterModels.parse_file(network_path)
            @test haskey(network["pump"]["1"], "efficiency_curve")

            # Add PUMP_LINEAR_POWER pump attributes here, even if they aren't used.
            map(x -> x["power_fixed"] = 0.0, values(network["pump"]))

            # Scale the expected power value and store as `expected_val`.
            transform_flow = WaterModels._calc_flow_per_unit_transform(network)
            transform_mass = WaterModels._calc_mass_per_unit_transform(network)
            transform_time = WaterModels._calc_time_per_unit_transform(network)
            transform_length = WaterModels._calc_length_per_unit_transform(network)
            expected_val = transform_mass(1.0) * transform_length(1.0)^2 /
                transform_time(1.0)^3 / transform_flow(1.0) * 100.0

            # Test that the power per unit flow is the expected value.
            map(x -> x["power_per_unit_flow"] = expected_val, values(network["pump"]))

            # Change the head curve form and recompute network bounds.
            map(x -> x["head_curve_form"] = head_curve_form, values(network["pump"]))
            WaterModels.recompute_bounds!(network)

            # Set up and solve an optimization problem with the pump data.
            set_flow_partitions!(network, 1.0, 1.0e-4)
            wm = instantiate_model(network, LRDWaterModel, build_wf)
            result = WaterModels.optimize_model!(wm, optimizer = _choose_solver(wm, ipopt, cbc))
            @test _is_valid_status(result["termination_status"])
        end
    end
end