@testset "src/core/pump.jl" begin
    for head_curve_form in [QUADRATIC, BEST_EFFICIENCY_POINT, EPANET]
        @testset "pump with efficiency curve: $(head_curve_form)" begin
            # Read in the initial network data and ensure an efficiency curve exists.
            network_path = "../test/data/epanet/snapshot/pump-hw-lps-eff-curve.inp"
            network = WaterModels.parse_file(network_path)
            @test haskey(network["pump"]["1"], "efficiency_curve")

            # Change the head curve form and recompute network bounds.
            map(x -> x["head_curve_form"] = head_curve_form, values(network["pump"]))
            WaterModels.recompute_bounds!(network)

            # Set up and solve an optimization problem with the pump data.
            ext = Dict(:pipe_breakpoints => 3, :pump_breakpoints => 3)
            wm = instantiate_model(network, LRDWaterModel, build_wf; ext = ext)
            result = WaterModels.optimize_model!(wm, optimizer = _make_juniper(wm, ipopt))
            @test result["termination_status"] == LOCALLY_SOLVED
        end
    end
end