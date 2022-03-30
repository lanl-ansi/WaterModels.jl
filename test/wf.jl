# Iterate over all possible WaterModels formulations.
for formulation in [NCWaterModel, NCDWaterModel, CRDWaterModel, LAWaterModel, LRDWaterModel, PWLRDWaterModel]
    @testset "Water Flow Problems (Single Network): $(formulation)" begin
        # Pipe, check valve, and shutoff valve tests will have the same physical solutions.
        for (component, name) in Dict("Pipe" => "pipe", "Check Valve" => "check_valve", "Shutoff Valve" => "shutoff_valve")
            @testset "Hazen-Williams Head Loss ($(component)): $(formulation)" begin
                network = WaterModels.parse_file("../test/data/epanet/snapshot/$(name)-hw-lps.inp")
                t_h = WaterModels._calc_head_per_unit_transform(network)
                t_q = WaterModels._calc_flow_per_unit_transform(network)
                set_flow_partitions_si!(network, 10.0, 1.0e-4)

                wm = instantiate_model(network, formulation, build_wf)
                result = WaterModels.optimize_model!(wm, optimizer = _choose_solver(wm, nlp_solver, milp_solver))

                @test _is_valid_status(result["termination_status"])
                @test isapprox(result["solution"]["node"]["1"]["h"], t_h(10.0), rtol = 1.0e-4)
                @test isapprox(result["solution"]["node"]["2"]["p"], t_h(7.89), rtol = 1.0e-2)
                @test isapprox(result["solution"]["pipe"]["1"]["q"], t_q(1.0), rtol = 1.0e-4)
            end
        end

        @testset "Hazen-Williams Head Loss (Negative Demand): $(formulation)" begin
            network = WaterModels.parse_file("../test/data/epanet/snapshot/negative_demand-hw-lps.inp")
            t_h = WaterModels._calc_head_per_unit_transform(network)
            t_q = WaterModels._calc_flow_per_unit_transform(network)
            set_flow_partitions_si!(network, 10.0, 1.0e-4)

            wm = instantiate_model(network, formulation, build_wf)
            result = WaterModels.optimize_model!(wm, optimizer = _choose_solver(wm, nlp_solver, milp_solver))

            @test _is_valid_status(result["termination_status"])
            @test isapprox(result["solution"]["node"]["1"]["h"], t_h(10.0), rtol = 1.0e-4)
            @test isapprox(result["solution"]["node"]["2"]["p"], t_h(8.27), rtol = 1.0e-2)
            @test isapprox(result["solution"]["node"]["3"]["p"], t_h(8.29), rtol = 1.0e-1)
            @test isapprox(result["solution"]["pipe"]["1"]["q"], t_q(0.9), rtol = 1.0e-4)
        end

        @testset "Head Loss (Regulator): $(formulation)" begin
            network = WaterModels.parse_file("../test/data/epanet/snapshot/prv-hw-lps.inp")
            t_h = WaterModels._calc_head_per_unit_transform(network)
            t_q = WaterModels._calc_flow_per_unit_transform(network)
            set_flow_partitions_si!(network, 10.0, 1.0e-4)

            wm = instantiate_model(network, formulation, build_wf)
            result = WaterModels.optimize_model!(wm, optimizer = _choose_solver(wm, nlp_solver, milp_solver))

            @test _is_valid_status(result["termination_status"])
            @test isapprox(result["solution"]["node"]["1"]["h"], t_h(10.0), rtol = 1.0e-4)
            @test isapprox(result["solution"]["node"]["2"]["p"], t_h(2.38), rtol = 1.0e-1)
            @test isapprox(result["solution"]["node"]["3"]["h"], t_h(2.00), rtol = 1.0e-2)
            @test isapprox(result["solution"]["node"]["3"]["p"], t_h(1.00), rtol = 1.0e-4)
            @test isapprox(result["solution"]["pipe"]["1"]["q"], t_q(2.0), rtol = 1.0e-4)
        end

        @testset "Darcy-Weisbach Head Loss: $(formulation)" begin
            network = WaterModels.parse_file("../test/data/epanet/snapshot/pipe-dw-lps.inp")
            t_h = WaterModels._calc_head_per_unit_transform(network)
            t_q = WaterModels._calc_flow_per_unit_transform(network)
            set_flow_partitions_si!(network, 10.0, 1.0e-4)

            wm = instantiate_model(network, formulation, build_wf)
            result = WaterModels.optimize_model!(wm, optimizer = _choose_solver(wm, nlp_solver, milp_solver))

            @test _is_valid_status(result["termination_status"])
            @test isapprox(result["solution"]["node"]["1"]["h"], t_h(10.0), rtol = 1.0e-4)
            @test isapprox(result["solution"]["node"]["2"]["p"], t_h(9.07), rtol = 1.0e-1)
            @test isapprox(result["solution"]["pipe"]["1"]["q"], t_q(1.0), rtol = 1.0e-4)
        end

        @testset "Short Pipe Dynamics: $(formulation)" begin
            network = WaterModels.parse_file("../test/data/epanet/snapshot/short-pipe-lps.inp")
            WaterModels.convert_short_pipes!(network)
            t_h = WaterModels._calc_head_per_unit_transform(network)
            t_q = WaterModels._calc_flow_per_unit_transform(network)
            set_flow_partitions_si!(network, 10.0, 1.0e-4)

            wm = instantiate_model(network, formulation, build_wf)
            result = WaterModels.optimize_model!(wm, optimizer = _choose_solver(wm, nlp_solver, milp_solver))

            @test _is_valid_status(result["termination_status"])
            @test isapprox(result["solution"]["node"]["1"]["h"], t_h(10.0), rtol = 1.0e-4)
            @test isapprox(result["solution"]["node"]["2"]["h"], t_h(10.0), rtol = 1.0e-4)
            @test isapprox(result["solution"]["short_pipe"]["1"]["q"], t_q(1.0), rtol = 1.0e-4)
        end

        @testset "Head Gain (Pump): $(formulation)" begin
            network = WaterModels.parse_file("../test/data/epanet/snapshot/pump-hw-lps.inp")
            t_h = WaterModels._calc_head_per_unit_transform(network)
            t_q = WaterModels._calc_flow_per_unit_transform(network)
            set_flow_partitions_si!(network, 10.0, 1.0e-4)

            wm = instantiate_model(network, formulation, build_wf)
            result = WaterModels.optimize_model!(wm, optimizer = _choose_solver(wm, nlp_solver, milp_solver))

            @test _is_valid_status(result["termination_status"])
            @test isapprox(result["solution"]["node"]["1"]["h"], t_h(10.0), rtol = 1.0e-3)
            @test isapprox(result["solution"]["pump"]["1"]["status"], 1.0, atol = 1.0e-3)
            @test result["solution"]["node"]["2"]["h"] > t_h(10.0)
        end

        @testset "Head Gain (Pump), Fix Pump Flow to Zero: $(formulation)" begin
            network = WaterModels.parse_file("../test/data/epanet/snapshot/pump-hw-lps.inp")
            network["pump"]["1"]["flow_min_forward"] = 0.0
            network["pump"]["1"]["flow_min"], network["pump"]["1"]["flow_max"] = 0.0, 0.0
            t_h = WaterModels._calc_head_per_unit_transform(network)
            t_q = WaterModels._calc_flow_per_unit_transform(network)
            set_flow_partitions_si!(network, 10.0, 1.0e-4)

            wm = instantiate_model(network, formulation, build_wf)
            result = WaterModels.optimize_model!(wm, optimizer = _choose_solver(wm, nlp_solver, milp_solver))
            @test !_is_valid_status(result["termination_status"])
        end

        @testset "Hazen-Williams Head Loss (Tank): $(formulation)" begin
            network = parse_file("../test/data/epanet/snapshot/tank-hw-lps.inp")
            t_h = WaterModels._calc_head_per_unit_transform(network)
            t_q = WaterModels._calc_flow_per_unit_transform(network)
            set_flow_partitions_si!(network, 10.0, 1.0e-4)

            wm = instantiate_model(network, formulation, build_wf)
            result = WaterModels.optimize_model!(wm, optimizer = _choose_solver(wm, nlp_solver, milp_solver))

            @test _is_valid_status(result["termination_status"])
            @test isapprox(result["solution"]["node"]["1"]["h"], t_h(20.0), rtol = 1.0e-3)
            @test isapprox(result["solution"]["node"]["2"]["h"], t_h(19.97), rtol = 2.5e-1)
            @test isapprox(result["solution"]["pipe"]["1"]["q"], t_q(0.1), rtol = 1.0e-3)
        end
    end

    @testset "Water Flow Problems (Multinetwork): $(formulation)" begin
        # Pipe, check valve, and shutoff valve tests will have the same physical solutions.
        for (component, name) in Dict("Pipe" => "pipe", "Check Valve" => "check_valve", "Shutoff Valve" => "shutoff_valve")
            @testset "Hazen-Williams Head Loss ($(component)): $(formulation)" begin
                network = WaterModels.parse_file("../test/data/epanet/multinetwork/$(name)-hw-lps.inp")
                network = WaterModels.make_multinetwork(network)
                t_h = WaterModels._calc_head_per_unit_transform(network)
                t_q = WaterModels._calc_flow_per_unit_transform(network)
                set_flow_partitions_si!(network, 10.0, 1.0e-4)

                wm = instantiate_model(network, formulation, build_mn_wf)
                result = WaterModels.optimize_model!(wm, optimizer = _choose_solver(wm, nlp_solver, milp_solver))

                @test _is_valid_status(result["termination_status"])
                @test isapprox(result["solution"]["nw"]["2"]["node"]["1"]["h"], t_h(10.0), rtol=1.0e-3)
                @test isapprox(result["solution"]["nw"]["1"]["node"]["2"]["p"], t_h(7.89), rtol=1.0e-1)
                @test isapprox(result["solution"]["nw"]["2"]["node"]["2"]["h"], t_h(9.42), rtol=1.0e-1)
                @test isapprox(result["solution"]["nw"]["3"]["node"]["2"]["h"], t_h(9.84), rtol=1.0e-1)
                @test isapprox(result["solution"]["nw"]["1"]["pipe"]["1"]["q"], t_q(1.0), rtol=1.0e-3)
                @test isapprox(result["solution"]["nw"]["2"]["pipe"]["1"]["q"], t_q(0.5), rtol=1.0e-3)
                @test isapprox(result["solution"]["nw"]["3"]["pipe"]["1"]["q"], t_q(0.25), rtol=1.0e-3)
            end
        end

        @testset "Hazen-Williams Head Loss (Negative Demand): $(formulation)" begin
            network = WaterModels.parse_file("../test/data/epanet/multinetwork/negative_demand-hw-lps.inp")
            network = WaterModels.make_multinetwork(network)
            t_h = WaterModels._calc_head_per_unit_transform(network)
            t_q = WaterModels._calc_flow_per_unit_transform(network)
            set_flow_partitions_si!(network, 10.0, 1.0e-4)

            wm = instantiate_model(network, formulation, build_mn_wf)
            result = WaterModels.optimize_model!(wm, optimizer = _choose_solver(wm, nlp_solver, milp_solver))

            @test _is_valid_status(result["termination_status"])
            @test isapprox(result["solution"]["nw"]["1"]["node"]["1"]["h"], t_h(10.0), rtol = 1.0e-3)
            @test isapprox(result["solution"]["nw"]["2"]["node"]["2"]["p"], t_h(8.27), rtol = 1.0e-1)
            @test isapprox(result["solution"]["nw"]["3"]["node"]["3"]["p"], t_h(8.29), rtol = 1.0e-1)
            @test isapprox(result["solution"]["nw"]["1"]["pipe"]["1"]["q"], t_q(0.9), rtol = 1.0e-3)
        end

        @testset "Head Loss (Regulator): $(formulation)" begin
            network = WaterModels.parse_file("../test/data/epanet/multinetwork/prv-hw-lps.inp")
            network = WaterModels.make_multinetwork(network)
            t_h = WaterModels._calc_head_per_unit_transform(network)
            t_q = WaterModels._calc_flow_per_unit_transform(network)
            set_flow_partitions_si!(network, 10.0, 1.0e-4)

            wm = instantiate_model(network, formulation, build_mn_wf)
            result = WaterModels.optimize_model!(wm, optimizer = _choose_solver(wm, nlp_solver, milp_solver))

            @test _is_valid_status(result["termination_status"])
            @test isapprox(result["solution"]["nw"]["1"]["node"]["2"]["h"], t_h(2.38), rtol = 2.5e-1)
            @test isapprox(result["solution"]["nw"]["2"]["node"]["2"]["h"], t_h(7.89), rtol = 2.5e-1)
            @test isapprox(result["solution"]["nw"]["3"]["node"]["2"]["h"], t_h(9.42), rtol = 2.5e-1)
            @test isapprox(result["solution"]["nw"]["2"]["node"]["3"]["h"], t_h(2.00), rtol = 1.0e-3)
            @test isapprox(result["solution"]["nw"]["2"]["node"]["1"]["h"], t_h(10.0), rtol = 1.0e-3)
            @test isapprox(result["solution"]["nw"]["1"]["pipe"]["1"]["q"], t_q(2.00), rtol = 1.0e-3)
            @test isapprox(result["solution"]["nw"]["3"]["pipe"]["1"]["q"], t_q(0.50), rtol = 1.0e-3)
            @test isapprox(result["solution"]["nw"]["1"]["regulator"]["2"]["q"], t_q(1.00), rtol = 1.0e-3)
            @test isapprox(result["solution"]["nw"]["3"]["regulator"]["2"]["q"], t_q(0.25), rtol = 1.0e-3)
        end

        @testset "Short Pipe Dynamics: $(formulation)" begin
            network = WaterModels.parse_file("../test/data/epanet/multinetwork/short-pipe-lps.inp")
            WaterModels.convert_short_pipes!(network)
            network = WaterModels.make_multinetwork(network)

            t_h = WaterModels._calc_head_per_unit_transform(network)
            t_q = WaterModels._calc_flow_per_unit_transform(network)
            set_flow_partitions_si!(network, 10.0, 1.0e-4)

            wm = instantiate_model(network, formulation, build_mn_wf)
            result = WaterModels.optimize_model!(wm, optimizer = _choose_solver(wm, nlp_solver, milp_solver))

            @test _is_valid_status(result["termination_status"])
            @test isapprox(result["solution"]["nw"]["1"]["node"]["2"]["h"], t_h(10.0), rtol = 1.0e-3)
            @test isapprox(result["solution"]["nw"]["1"]["short_pipe"]["1"]["q"], t_q(1.0), rtol = 1.0e-3)
            @test isapprox(result["solution"]["nw"]["3"]["node"]["2"]["h"], t_h(10.0), rtol = 1.0e-3)
            @test isapprox(result["solution"]["nw"]["3"]["short_pipe"]["1"]["q"], t_q(0.25), rtol = 1.0e-3)
        end

        @testset "Head Gain (Pump): $(formulation)" begin
            network = WaterModels.parse_file("../test/data/epanet/multinetwork/pump-hw-lps.inp")
            network = WaterModels.make_multinetwork(network)
            t_h = WaterModels._calc_head_per_unit_transform(network)
            t_q = WaterModels._calc_flow_per_unit_transform(network)
            set_flow_partitions_si!(network, 10.0, 1.0e-4)

            wm = instantiate_model(network, formulation, build_mn_wf)
            result = WaterModels.optimize_model!(wm, optimizer = _choose_solver(wm, nlp_solver, milp_solver))

            @test _is_valid_status(result["termination_status"])
            @test isapprox(result["solution"]["nw"]["1"]["pump"]["1"]["q"], t_q(0.125), rtol = 1.0e-3)
            @test isapprox(result["solution"]["nw"]["3"]["pump"]["1"]["q"], t_q(0.03125), rtol = 1.0e-3)
            @test isapprox(result["solution"]["nw"]["1"]["pump"]["1"]["status"], 1.0, atol = 1.0e-3)
            @test result["solution"]["nw"]["1"]["pump"]["1"]["g"] > 0.0
            @test result["solution"]["nw"]["3"]["pump"]["1"]["g"] > 0.0
        end

        @testset "Hazen-Williams Head Loss (Tank): $(formulation)" begin
            network = WaterModels.parse_file("../test/data/epanet/multinetwork/tank-hw-lps.inp")
            network = WaterModels.make_multinetwork(network)
            t_h = WaterModels._calc_head_per_unit_transform(network)
            t_q = WaterModels._calc_flow_per_unit_transform(network)
            set_flow_partitions_si!(network, 10.0, 1.0e-4)

            wm = instantiate_model(network, formulation, build_mn_wf)
            result = WaterModels.optimize_model!(wm, optimizer = _choose_solver(wm, nlp_solver, milp_solver))

            @test _is_valid_status(result["termination_status"])
            @test isapprox(result["solution"]["nw"]["1"]["node"]["1"]["h"], t_h(60.00), rtol = 1.0e-3)
            @test isapprox(result["solution"]["nw"]["1"]["node"]["2"]["h"], t_h(59.42), rtol = 1.0e-1)
            @test isapprox(result["solution"]["nw"]["3"]["node"]["1"]["h"], t_h(25.62), rtol = 1.0e-3)
            @test isapprox(result["solution"]["nw"]["3"]["node"]["2"]["h"], t_h(25.58), rtol = 1.0e-1)
            @test isapprox(result["solution"]["nw"]["1"]["pipe"]["1"]["q"], t_q(0.500), rtol = 1.0e-3)
            @test isapprox(result["solution"]["nw"]["3"]["pipe"]["1"]["q"], t_q(0.125), rtol = 1.0e-3)
        end
    end
end


@testset "solve_wf" begin
    network = WaterModels.parse_file("../test/data/epanet/snapshot/pipe-hw-lps.inp")
    set_flow_partitions_si!(network, 10.0, 1.0e-4)
    
    result = WaterModels.solve_wf(network, LRDWaterModel, milp_solver)
    result = WaterModels.run_wf(network, LRDWaterModel, milp_solver)
    @test _is_valid_status(result["termination_status"])
end


@testset "solve_mn_wf" begin
    @testset "multinetwork data" begin
        network = WaterModels.parse_file("../test/data/epanet/multinetwork/pipe-hw-lps.inp")
        network_mn = WaterModels.make_multinetwork(network)
        set_flow_partitions_si!(network_mn, 10.0, 1.0e-4)

        result = WaterModels.solve_mn_wf(network_mn, LRDWaterModel, milp_solver)
        result = WaterModels.run_mn_wf(network_mn, LRDWaterModel, milp_solver)
        @test _is_valid_status(result["termination_status"])
    end

    @testset "single network data" begin
        network = WaterModels.parse_file("../test/data/epanet/snapshot/pipe-hw-lps.inp")
        network_mn = WaterModels.make_multinetwork(network)
        set_flow_partitions_si!(network_mn, 10.0, 1.0e-4)

        result = WaterModels.solve_mn_wf(network_mn, LRDWaterModel, milp_solver)
        result = WaterModels.run_mn_wf(network_mn, LRDWaterModel, milp_solver)
        @test _is_valid_status(result["termination_status"])
    end
end


@testset "solve_mn_wf_switching" begin
    network = WaterModels.parse_file("../test/data/epanet/multinetwork/pump-hw-lps.inp")
    network_mn = WaterModels.make_multinetwork(network)
    set_flow_partitions_si!(network_mn, 10.0, 1.0e-4)

    result = WaterModels.solve_mn_wf_switching(network_mn, LRDWaterModel, milp_solver)
    result = WaterModels.run_mn_wf_switching(network_mn, LRDWaterModel, milp_solver)
    @test _is_valid_status(result["termination_status"])
end


@testset "solve_wf (with symmetric pumps)" begin
    network = WaterModels.parse_file("../test/data/epanet/snapshot/pump-hw-lps-multiple.inp"; skip_correct = true)
    modifications = WaterModels.parse_file("../test/data/json/pump-multiple-group.json"; skip_correct = true)
    _IM.update_data!(network, modifications)
    WaterModels.correct_network_data!(network)

    set_flow_partitions_si!(network, 10.0, 1.0e-4)
    result = WaterModels.solve_wf(network, LRDWaterModel, milp_solver)
    @test _is_valid_status(result["termination_status"])
end


@testset "solve_mn_wf (with symmetric pumps)" begin
    network = parse_file("../test/data/epanet/multinetwork/pump-hw-lps-multiple.inp"; skip_correct = true)
    modifications = WaterModels.parse_file("../test/data/json/pump-multiple-group.json"; skip_correct = true)
    _IM.update_data!(network, modifications)
    WaterModels.correct_network_data!(network)
    network_mn = WaterModels.make_multinetwork(network)
    set_flow_partitions_si!(network_mn, 10.0, 1.0e-4)

    result = WaterModels.solve_mn_wf(network_mn, LRDWaterModel, milp_solver)
    @test _is_valid_status(result["termination_status"])
end