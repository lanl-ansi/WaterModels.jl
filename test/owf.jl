# Iterate over all possible WaterModels formulations.
for formulation in [NCWaterModel, NCDWaterModel, CRDWaterModel, LAWaterModel, LRDWaterModel, PWLRDWaterModel]
    # Set a generic extensions dictionary for testing purposes.
    ext = Dict(:pipe_breakpoints => 3, :pump_breakpoints => 3)

    @testset "Optimal Water Flow Problems (Single Network): $(formulation)" begin
        network = WaterModels.parse_file("../test/data/epanet/snapshot/pump-hw-lps.inp")
        wm = instantiate_model(network, formulation, build_owf; ext = ext)
        result = WaterModels.optimize_model!(wm, optimizer = _choose_solver(wm, ipopt, cbc))

        @test _is_valid_status(result["termination_status"])
        @test isapprox(result["solution"]["node"]["1"]["h"], 10.0, rtol = 1.0e-3)
        @test isapprox(result["solution"]["pump"]["1"]["status"], 1.0, atol = 1.0e-3)
        @test result["solution"]["node"]["2"]["h"] > 10.0

        # In some relaxations, the objective may be near zero.
        @test result["objective"] >= -1.0e-6
    end

    @testset "Optimal Water Flow Problems (Single Network, Best Efficiency Point): $(formulation)" begin
        network = WaterModels.parse_file("../test/data/epanet/snapshot/pump-hw-lps.inp")
        map(x -> x["head_curve_form"] = WaterModels.BEST_EFFICIENCY_POINT, values(network["pump"]))
        WaterModels.recompute_bounds!(network) # Recompute component bounds after the above changes.
        wm = instantiate_model(network, formulation, build_owf; ext = ext)
        result = WaterModels.optimize_model!(wm, optimizer = _choose_solver(wm, ipopt, cbc))

        @test _is_valid_status(result["termination_status"])
        @test isapprox(result["solution"]["node"]["1"]["h"], 10.0, rtol = 1.0e-3)
        @test isapprox(result["solution"]["pump"]["1"]["status"], 1.0, atol = 1.0e-3)
        @test result["solution"]["node"]["2"]["h"] > 10.0

        # In some relaxations, the objective may be near zero.
        @test result["objective"] >= -1.0e-6
    end

    @testset "Optimal Water Flow Problems (Single Network, EPANET Pump Formulation): $(formulation)" begin
        network = WaterModels.parse_file("../test/data/epanet/snapshot/pump-hw-lps.inp")
        map(x -> x["head_curve_form"] = WaterModels.EPANET, values(network["pump"]))
        WaterModels.recompute_bounds!(network) # Recompute component bounds after the above changes.

        if !(formulation <: AbstractNonlinearModel)
            wm = instantiate_model(network, formulation, build_owf; ext = ext)
            result = WaterModels.optimize_model!(wm, optimizer = _choose_solver(wm, ipopt, cbc))

            @test _is_valid_status(result["termination_status"])
            @test isapprox(result["solution"]["node"]["1"]["h"], 10.0, rtol = 1.0e-3)
            @test isapprox(result["solution"]["pump"]["1"]["status"], 1.0, atol = 1.0e-3)
            @test result["solution"]["node"]["2"]["h"] > 10.0

            # In some relaxations, the objective may be near zero.
            @test result["objective"] >= -1.0e-6
        end
    end

    @testset "Optimal Water Flow Problems (Multinetwork): $(formulation)" begin
        network = WaterModels.parse_file("../test/data/epanet/multinetwork/owf-hw-lps.inp")
        network_mn = WaterModels.make_multinetwork(network)
        wm = instantiate_model(network_mn, formulation, build_mn_owf; ext = ext)
        result = WaterModels.optimize_model!(wm, optimizer = _choose_solver(wm, ipopt, cbc))

        @test _is_valid_status(result["termination_status"])
    end
end


@testset "run_owf" begin
    network = WaterModels.parse_file("../test/data/epanet/snapshot/pump-hw-lps.inp")
    result = WaterModels.run_owf(network, LRDWaterModel, cbc)
    @test result["termination_status"] == OPTIMAL
end


@testset "run_mn_owf" begin
    network = WaterModels.parse_file("../test/data/epanet/multinetwork/owf-hw-lps.inp")
    network_mn = WaterModels.make_multinetwork(network)
    result = WaterModels.run_mn_owf(network_mn, LRDWaterModel, cbc)
    @test result["termination_status"] == OPTIMAL
end