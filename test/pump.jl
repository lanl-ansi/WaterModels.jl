@testset "src/core/pump.jl" begin
    for head_curve_form in [PUMP_QUADRATIC, PUMP_BEST_EFFICIENCY_POINT, PUMP_EPANET, PUMP_LINEAR_POWER]
        @testset "pump with efficiency curve: $(head_curve_form)" begin
            # Read in the initial network data and ensure an efficiency curve exists.
            network_path = "../test/data/epanet/snapshot/pump-hw-lps-eff-curve.inp"
            network = WaterModels.parse_file(network_path)
            @test haskey(network["pump"]["1"], "efficiency_curve")

            # Add PUMP_LINEAR_POWER pump attributes here, even if they aren't used.
            map(x -> x["power_fixed"] = 0.0, values(network["pump"]))
            map(x -> x["power_per_unit_flow"] = 100.0, values(network["pump"]))

            # Change the head curve form and recompute network bounds.
            map(x -> x["head_curve_form"] = head_curve_form, values(network["pump"]))
            WaterModels.recompute_bounds!(network)

            # Set up and solve an optimization problem with the pump data.
            ext = Dict(:pipe_breakpoints => 3, :pump_breakpoints => 3)
            wm = instantiate_model(network, LRDWaterModel, build_wf; ext = ext)
            result = WaterModels.optimize_model!(wm, optimizer = _choose_solver(wm, ipopt, cbc))
            @test _is_valid_status(result["termination_status"])
        end
    end
end