# Iterate over all possible WaterModels formulations.
for formulation in [NCWaterModel, NCDWaterModel, CRDWaterModel, LAWaterModel, LRDWaterModel, PWLRDWaterModel]
    @testset "Optimal Water Flow Problems (Single Network): $(formulation)" begin
        network = WaterModels.parse_file("../test/data/epanet/snapshot/pump-hw-lps.inp")
        set_flow_partitions_si!(network, 10.0, 1.0e-4)
        t_h = WaterModels._calc_head_per_unit_transform(network)

        wm = instantiate_model(network, formulation, build_owf)
        result = WaterModels.optimize_model!(wm, optimizer = _choose_solver(wm, nlp_solver, milp_solver))

        @test _is_valid_status(result["termination_status"])
        @test isapprox(result["solution"]["node"]["1"]["h"], t_h(10.0), rtol = 1.0e-3)
        @test isapprox(result["solution"]["pump"]["1"]["status"], 1.0, atol = 1.0e-3)
        @test result["solution"]["node"]["2"]["h"] > t_h(10.0)

        # In some relaxations, the objective may be near zero.
        @test result["objective"] >= -1.0e-6
    end

    @testset "Optimal Water Flow Problems (Single Network, Best Efficiency Point): $(formulation)" begin
        network = WaterModels.parse_file("../test/data/epanet/snapshot/pump-hw-lps.inp")
        t_h = WaterModels._calc_head_per_unit_transform(network)

        map(x -> x["pump_type"] = WaterModels.PUMP_BEST_EFFICIENCY_POINT, values(network["pump"]))
        WaterModels.recompute_bounds!(network) # Recompute component bounds after the above changes.
        set_flow_partitions_si!(network, 10.0, 1.0e-4)

        wm = instantiate_model(network, formulation, build_owf)
        result = WaterModels.optimize_model!(wm, optimizer = _choose_solver(wm, nlp_solver, milp_solver))

        @test _is_valid_status(result["termination_status"])
        @test isapprox(result["solution"]["node"]["1"]["h"], t_h(10.0), rtol = 1.0e-3)
        @test isapprox(result["solution"]["pump"]["1"]["status"], 1.0, atol = 1.0e-3)
        @test result["solution"]["node"]["2"]["h"] > t_h(10.0)

        # In some relaxations, the objective may be near zero.
        @test result["objective"] >= -1.0e-6
    end

    @testset "Optimal Water Flow Problems (Single Network, EPANET Pump Formulation): $(formulation)" begin
        network = WaterModels.parse_file("../test/data/epanet/snapshot/pump-hw-lps.inp")
        t_h = WaterModels._calc_head_per_unit_transform(network)

        map(x -> x["pump_type"] = PUMP_EPANET, values(network["pump"]))
        WaterModels.recompute_bounds!(network) # Recompute component bounds after the above changes.
        set_flow_partitions_si!(network, 10.0, 1.0e-4)

        if !(formulation <: AbstractNonlinearModel)
            wm = instantiate_model(network, formulation, build_owf)
            result = WaterModels.optimize_model!(wm, optimizer = _choose_solver(wm, nlp_solver, milp_solver))

            @test _is_valid_status(result["termination_status"])
            @test isapprox(result["solution"]["node"]["1"]["h"], t_h(10.0), rtol = 1.0e-3)
            @test isapprox(result["solution"]["pump"]["1"]["status"], 1.0, atol = 1.0e-3)
            @test result["solution"]["node"]["2"]["h"] > t_h(10.0)

            # In some relaxations, the objective may be near zero.
            @test result["objective"] >= -1.0e-6
        end
    end

    @testset "Optimal Water Flow Problems (Multinetwork): $(formulation)" begin
        network = WaterModels.parse_file("../test/data/epanet/multinetwork/owf-hw-lps.inp")
        network_mn = WaterModels.make_multinetwork(network)
        set_flow_partitions_si!(network_mn, 10.0, 1.0e-4)

        wm = instantiate_model(network_mn, formulation, build_mn_owf)
        result = WaterModels.optimize_model!(wm, optimizer = _choose_solver(wm, nlp_solver, milp_solver))

        @test _is_valid_status(result["termination_status"])
    end
end


@testset "solve_owf" begin
    network = WaterModels.parse_file("../test/data/epanet/snapshot/pump-hw-lps.inp")
    set_flow_partitions_si!(network, 10.0, 1.0e-4)

    result = WaterModels.solve_owf(network, LRDWaterModel, milp_solver)
    result = WaterModels.run_owf(network, LRDWaterModel, milp_solver)
    @test result["termination_status"] == OPTIMAL
end


@testset "solve_mn_owf" begin
    network = WaterModels.parse_file("../test/data/epanet/multinetwork/owf-hw-lps.inp")
    network_mn = WaterModels.make_multinetwork(network)
    set_flow_partitions_si!(network_mn, 10.0, 1.0e-4)

    result = WaterModels.solve_mn_owf(network_mn, LRDWaterModel, milp_solver)
    result = WaterModels.run_mn_owf(network_mn, LRDWaterModel, milp_solver)
    @test result["termination_status"] == OPTIMAL
end


@testset "solve_mn_owf_switching" begin
    network = WaterModels.parse_file("../test/data/epanet/multinetwork/owf-hw-lps.inp")
    network_mn = WaterModels.make_multinetwork(network)
    set_flow_partitions_si!(network_mn, 10.0, 1.0e-4)

    result = WaterModels.solve_mn_owf_switching(network_mn, LRDWaterModel, milp_solver)
    result = WaterModels.run_mn_owf_switching(network_mn, LRDWaterModel, milp_solver)
    @test result["termination_status"] == OPTIMAL
end