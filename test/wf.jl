# Iterate over all possible WaterModels formulations.
for formulation in [NCWaterModel, NCDWaterModel, CRDWaterModel, LAWaterModel, LRDWaterModel, PWLRDWaterModel]
    # Set a generic extensions dictionary for testing purposes.
    ext = Dict(:pipe_breakpoints => 5, :pump_breakpoints => 5)

    @testset "Water Flow Problems (Single Network): $(formulation)" begin
        # Pipe, check valve, and shutoff valve tests will have the same physical solutions.
        for (component, name) in Dict("Pipe" => "pipe", "Check Valve" => "check_valve", "Shutoff Valve" => "shutoff_valve")
            @testset "Hazen-Williams Head Loss ($(component)): $(formulation)" begin
                network = WaterModels.parse_file("../test/data/epanet/snapshot/$(name)-hw-lps.inp")
                wm = instantiate_model(network, formulation, build_wf; ext = ext)
                result = WaterModels.optimize_model!(wm, optimizer = _choose_solver(wm, ipopt, cbc))

                @test _is_valid_status(result["termination_status"])
                @test isapprox(result["solution"]["node"]["1"]["h"], 10.0, rtol = 1.0e-4)
                @test isapprox(result["solution"]["node"]["2"]["p"], 7.89, rtol = 1.0e-2)
                @test isapprox(result["solution"]["pipe"]["1"]["q"], 1.0, rtol = 1.0e-4)
            end
        end

        @testset "Hazen-Williams Head Loss (Negative Demand): $(formulation)" begin
            network = WaterModels.parse_file("../test/data/epanet/snapshot/negative_demand-hw-lps.inp")
            wm = instantiate_model(network, formulation, build_wf; ext = ext)
            result = WaterModels.optimize_model!(wm, optimizer = _choose_solver(wm, ipopt, cbc))

            @test _is_valid_status(result["termination_status"])
            @test isapprox(result["solution"]["node"]["1"]["h"], 10.0, rtol = 1.0e-4)
            @test isapprox(result["solution"]["node"]["2"]["p"], 8.27, rtol = 1.0e-2)
            @test isapprox(result["solution"]["node"]["3"]["p"], 8.29, rtol = 1.0e-1)
            @test isapprox(result["solution"]["pipe"]["1"]["q"], 0.9, rtol = 1.0e-4)
        end

        @testset "Head Loss (Regulator): $(formulation)" begin
            network = WaterModels.parse_file("../test/data/epanet/snapshot/prv-hw-lps.inp")
            wm = instantiate_model(network, formulation, build_wf; ext = ext)
            result = WaterModels.optimize_model!(wm, optimizer = _choose_solver(wm, ipopt, cbc))

            @test _is_valid_status(result["termination_status"])
            @test isapprox(result["solution"]["node"]["1"]["h"], 10.0, rtol = 1.0e-4)
            @test isapprox(result["solution"]["node"]["2"]["p"], 2.38, rtol = 1.0e-1)
            @test isapprox(result["solution"]["node"]["3"]["h"], 2.00, rtol = 1.0e-2)
            @test isapprox(result["solution"]["node"]["3"]["p"], 1.00, rtol = 1.0e-4)
            @test isapprox(result["solution"]["pipe"]["1"]["q"], 2.0, rtol = 1.0e-4)
        end

        @testset "Short Pipe Dynamics: $(formulation)" begin
            network = WaterModels.parse_file("../test/data/epanet/snapshot/short-pipe-lps.inp")
            wm = instantiate_model(network, formulation, build_wf; ext = ext)
            result = WaterModels.optimize_model!(wm, optimizer = _choose_solver(wm, ipopt, cbc))

            @test _is_valid_status(result["termination_status"])
            @test isapprox(result["solution"]["node"]["1"]["h"], 10.0, rtol = 1.0e-4)
            @test isapprox(result["solution"]["node"]["2"]["h"], 10.0, rtol = 1.0e-4)
            @test isapprox(result["solution"]["short_pipe"]["1"]["q"], 1.0, rtol = 1.0e-4)
        end

        @testset "Head Gain (Pump): $(formulation)" begin
            network = WaterModels.parse_file("../test/data/epanet/snapshot/pump-hw-lps.inp")
            wm = instantiate_model(network, formulation, build_wf; ext = ext)
            result = WaterModels.optimize_model!(wm, optimizer = _choose_solver(wm, ipopt, cbc))

            @test _is_valid_status(result["termination_status"])
            @test isapprox(result["solution"]["node"]["1"]["h"], 10.0, rtol = 1.0e-3)
            @test isapprox(result["solution"]["pump"]["1"]["status"], 1.0, atol = 1.0e-3)
            @test result["solution"]["node"]["2"]["h"] > 10.0
        end

        @testset "Hazen-Williams Head Loss (Tank): $(formulation)" begin
            network = parse_file("../test/data/epanet/snapshot/tank-hw-lps.inp")
            wm = instantiate_model(network, formulation, build_wf; ext = ext)
            result = WaterModels.optimize_model!(wm, optimizer = _choose_solver(wm, ipopt, cbc))

            @test _is_valid_status(result["termination_status"])
            @test isapprox(result["solution"]["node"]["1"]["h"], 20.0, rtol = 1.0e-3)
            @test isapprox(result["solution"]["node"]["2"]["h"], 17.89, rtol = 2.5e-1)
            @test isapprox(result["solution"]["pipe"]["1"]["q"], 1.0, rtol = 1.0e-3)
        end
    end

    @testset "Water Flow Problems (Multinetwork): $(formulation)" begin
        # Pipe, check valve, and shutoff valve tests will have the same physical solutions.
        for (component, name) in Dict("Pipe" => "pipe", "Check Valve" => "check_valve", "Shutoff Valve" => "shutoff_valve")
            @testset "Hazen-Williams Head Loss ($(component)): $(formulation)" begin
                network = WaterModels.parse_file("../test/data/epanet/multinetwork/$(name)-hw-lps.inp")
                network = WaterModels.make_multinetwork(network)
                wm = instantiate_model(network, formulation, build_mn_wf; ext = ext)
                result = WaterModels.optimize_model!(wm, optimizer = _choose_solver(wm, ipopt, cbc))

                @test _is_valid_status(result["termination_status"])
                @test isapprox(result["solution"]["nw"]["2"]["node"]["1"]["h"], 10.0, rtol=1.0e-3)
                @test isapprox(result["solution"]["nw"]["1"]["node"]["2"]["p"], 7.89, rtol=1.0e-1)
                @test isapprox(result["solution"]["nw"]["2"]["node"]["2"]["h"], 9.42, rtol=1.0e-1)
                @test isapprox(result["solution"]["nw"]["3"]["node"]["2"]["h"], 9.84, rtol=1.0e-1)
                @test isapprox(result["solution"]["nw"]["1"]["pipe"]["1"]["q"], 1.0, rtol=1.0e-3)
                @test isapprox(result["solution"]["nw"]["2"]["pipe"]["1"]["q"], 0.5, rtol=1.0e-3)
                @test isapprox(result["solution"]["nw"]["3"]["pipe"]["1"]["q"], 0.25, rtol=1.0e-3)
            end
        end

        @testset "Hazen-Williams Head Loss (Negative Demand): $(formulation)" begin
            network = WaterModels.parse_file("../test/data/epanet/multinetwork/negative_demand-hw-lps.inp")
            network = WaterModels.make_multinetwork(network)
            wm = instantiate_model(network, formulation, build_mn_wf; ext = ext)
            result = WaterModels.optimize_model!(wm, optimizer = _choose_solver(wm, ipopt, cbc))

            @test _is_valid_status(result["termination_status"])
            @test isapprox(result["solution"]["nw"]["1"]["node"]["1"]["h"], 10.0, rtol = 1.0e-3)
            @test isapprox(result["solution"]["nw"]["2"]["node"]["2"]["p"], 8.27, rtol = 1.0e-1)
            @test isapprox(result["solution"]["nw"]["3"]["node"]["3"]["p"], 8.29, rtol = 1.0e-1)
            @test isapprox(result["solution"]["nw"]["1"]["pipe"]["1"]["q"], 0.9, rtol = 1.0e-3)
        end

        @testset "Head Loss (Regulator): $(formulation)" begin
            network = WaterModels.parse_file("../test/data/epanet/multinetwork/prv-hw-lps.inp")
            network = WaterModels.make_multinetwork(network)
            wm = instantiate_model(network, formulation, build_mn_wf; ext = ext)
            result = WaterModels.optimize_model!(wm, optimizer = _choose_solver(wm, ipopt, cbc))

            @test _is_valid_status(result["termination_status"])
            @test isapprox(result["solution"]["nw"]["1"]["node"]["2"]["h"], 2.38, rtol = 2.5e-1)
            @test isapprox(result["solution"]["nw"]["2"]["node"]["2"]["h"], 7.89, rtol = 2.5e-1)
            @test isapprox(result["solution"]["nw"]["3"]["node"]["2"]["h"], 9.42, rtol = 2.5e-1)
            @test isapprox(result["solution"]["nw"]["2"]["node"]["3"]["h"], 2.00, rtol = 1.0e-3)
            @test isapprox(result["solution"]["nw"]["2"]["node"]["1"]["h"], 10.0, rtol = 1.0e-3)
            @test isapprox(result["solution"]["nw"]["1"]["pipe"]["1"]["q"], 2.00, rtol = 1.0e-3)
            @test isapprox(result["solution"]["nw"]["3"]["pipe"]["1"]["q"], 0.50, rtol = 1.0e-3)
            @test isapprox(result["solution"]["nw"]["1"]["regulator"]["2"]["q"], 1.00, rtol = 1.0e-3)
            @test isapprox(result["solution"]["nw"]["3"]["regulator"]["2"]["q"], 0.25, rtol = 1.0e-3)
        end

        @testset "Short Pipe Dynamics: $(formulation)" begin
            network = WaterModels.parse_file("../test/data/epanet/multinetwork/short-pipe-lps.inp")
            network = WaterModels.make_multinetwork(network)
            wm = instantiate_model(network, formulation, build_mn_wf; ext = ext)
            result = WaterModels.optimize_model!(wm, optimizer = _choose_solver(wm, ipopt, cbc))

            @test _is_valid_status(result["termination_status"])
            @test isapprox(result["solution"]["nw"]["1"]["node"]["2"]["h"], 10.0, rtol = 1.0e-3)
            @test isapprox(result["solution"]["nw"]["1"]["short_pipe"]["1"]["q"], 1.0, rtol = 1.0e-3)
            @test isapprox(result["solution"]["nw"]["3"]["node"]["2"]["h"], 10.0, rtol = 1.0e-3)
            @test isapprox(result["solution"]["nw"]["3"]["short_pipe"]["1"]["q"], 0.25, rtol = 1.0e-3)
        end

        @testset "Head Gain (Pump): $(formulation)" begin
            network = WaterModels.parse_file("../test/data/epanet/multinetwork/pump-hw-lps.inp")
            network = WaterModels.make_multinetwork(network)
            wm = instantiate_model(network, formulation, build_mn_wf; ext = ext)
            result = WaterModels.optimize_model!(wm, optimizer = _choose_solver(wm, ipopt, cbc))

            @test _is_valid_status(result["termination_status"])
            @test isapprox(result["solution"]["nw"]["1"]["pump"]["1"]["q"], 0.125, rtol = 1.0e-3)
            @test isapprox(result["solution"]["nw"]["3"]["pump"]["1"]["q"], 0.03125, rtol = 1.0e-3)
            @test isapprox(result["solution"]["nw"]["1"]["pump"]["1"]["status"], 1.0, atol = 1.0e-3)
            @test result["solution"]["nw"]["1"]["pump"]["1"]["g"] > 0.0
            @test result["solution"]["nw"]["3"]["pump"]["1"]["g"] > 0.0
        end

        @testset "Hazen-Williams Head Loss (Tank): $(formulation)" begin
            network = WaterModels.parse_file("../test/data/epanet/multinetwork/tank-hw-lps.inp")
            network = WaterModels.make_multinetwork(network)
            wm = instantiate_model(network, formulation, build_mn_wf; ext = ext)
            result = WaterModels.optimize_model!(wm, optimizer = _choose_solver(wm, ipopt, cbc))

            @test _is_valid_status(result["termination_status"])
            @test isapprox(result["solution"]["nw"]["1"]["node"]["1"]["h"], 60.00, rtol = 1.0e-3)
            @test isapprox(result["solution"]["nw"]["1"]["node"]["2"]["h"], 59.42, rtol = 1.0e-1)
            @test isapprox(result["solution"]["nw"]["3"]["node"]["1"]["h"], 25.62, rtol = 1.0e-3)
            @test isapprox(result["solution"]["nw"]["3"]["node"]["2"]["h"], 25.58, rtol = 1.0e-1)
            @test isapprox(result["solution"]["nw"]["1"]["pipe"]["1"]["q"], 0.500, rtol = 1.0e-3)
            @test isapprox(result["solution"]["nw"]["3"]["pipe"]["1"]["q"], 0.125, rtol = 1.0e-3)
        end
    end
end


@testset "run_wf" begin
    network = WaterModels.parse_file("../test/data/epanet/snapshot/pipe-hw-lps.inp")
    result = WaterModels.run_wf(network, LRDWaterModel, cbc)
    @test _is_valid_status(result["termination_status"])
end


@testset "run_mn_wf" begin
    network = WaterModels.parse_file("../test/data/epanet/multinetwork/pipe-hw-lps.inp")
    network_mn = WaterModels.make_multinetwork(network)
    result = WaterModels.run_mn_wf(network_mn, LRDWaterModel, cbc)
    @test _is_valid_status(result["termination_status"])
end