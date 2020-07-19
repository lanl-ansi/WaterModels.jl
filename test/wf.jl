@testset "Water Flow Problems (Single Network)" begin
    # Pipe, check valve, and shutoff valve tests will have the same physical solutions.
    for (component, name) in Dict("Pipe"=>"pipe", "Check Valve"=>"check_valve", "Shutoff Valve"=>"shutoff_valve")
        @testset "Hazen-Williams Head Loss ($(component))" begin
            network = WaterModels.parse_file("../test/data/epanet/snapshot/$(name)-hw-lps.inp")

            wm = instantiate_model(network, NLPWaterModel, build_wf)
            result = WaterModels.optimize_model!(wm, optimizer=_make_juniper(wm, ipopt))
            @test result["termination_status"] == LOCALLY_SOLVED
            @test isapprox(result["solution"]["node"]["1"]["h"], 10.0, rtol=1.0e-3)
            @test isapprox(result["solution"]["node"]["2"]["p"], 7.89, rtol=1.0e-3)
            @test isapprox(result["solution"][name]["1"]["q"], 1.0, rtol=1.0e-3)

            wm = instantiate_model(network, MICPRWaterModel, build_wf)
            result = WaterModels.optimize_model!(wm, optimizer=_make_juniper(wm, ipopt))
            @test result["termination_status"] == LOCALLY_SOLVED
            @test isapprox(result["solution"]["node"]["1"]["h"], 10.0, rtol=1.0e-3)
            @test result["solution"]["node"]["2"]["p"] <= 7.89
            @test isapprox(result["solution"][name]["1"]["q"], 1.0, rtol=1.0e-3)

            result = run_wf(network, MILPWaterModel, cbc, ext=Dict(:pipe_breakpoints=>3))
            @test result["termination_status"] == OPTIMAL
            @test isapprox(result["solution"]["node"]["1"]["h"], 10.0, rtol=1.0e-3)
            @test isapprox(result["solution"]["node"]["2"]["p"], 7.89, rtol=1.0e-3)
            @test isapprox(result["solution"][name]["1"]["q"], 1.0, rtol=1.0e-3)

            result = run_wf(network, MILPRWaterModel, cbc, ext=Dict(:pipe_breakpoints=>3))
            @test result["termination_status"] == OPTIMAL
            @test isapprox(result["solution"]["node"]["1"]["h"], 10.0, rtol=1.0e-3)
            @test isapprox(result["solution"][name]["1"]["q"], 1.0, rtol=1.0e-3)
        end
    end

    @testset "Hazen-Williams Head Loss (Negative Demand)" begin
        network = WaterModels.parse_file("../test/data/epanet/snapshot/negative_demand-hw-lps.inp")

        wm = instantiate_model(network, NLPWaterModel, build_wf)
        result = WaterModels.optimize_model!(wm, optimizer=_make_juniper(wm, ipopt))
        @test result["termination_status"] == LOCALLY_SOLVED
        @test isapprox(result["solution"]["node"]["1"]["h"], 10.0, rtol=1.0e-3)
        @test isapprox(result["solution"]["node"]["2"]["p"], 8.27, rtol=1.0e-3)
        @test isapprox(result["solution"]["node"]["3"]["p"], 8.29, rtol=1.0e-3)
        @test isapprox(result["solution"]["pipe"]["1"]["q"], 0.9, rtol=1.0e-3)

        wm = instantiate_model(network, MICPRWaterModel, build_wf)
        result = WaterModels.optimize_model!(wm, optimizer=_make_juniper(wm, ipopt))
        @test result["termination_status"] == LOCALLY_SOLVED
        @test isapprox(result["solution"]["node"]["1"]["h"], 10.0, rtol=1.0e-3)
        @test isapprox(result["solution"]["pipe"]["1"]["q"], 0.9, rtol=1.0e-3)

        result = run_wf(network, MILPWaterModel, cbc, ext=Dict(:pipe_breakpoints=>4))
        @test result["termination_status"] == OPTIMAL
        @test isapprox(result["solution"]["node"]["1"]["h"], 10.0, rtol=1.0e-3)
        @test isapprox(result["solution"]["node"]["2"]["p"], 8.27, rtol=1.0e-2)
        @test isapprox(result["solution"]["node"]["3"]["p"], 8.29, rtol=1.0e-3)
        @test isapprox(result["solution"]["pipe"]["1"]["q"], 0.9, rtol=1.0e-3)

        result = run_wf(network, MILPRWaterModel, cbc, ext=Dict(:pipe_breakpoints=>3))
        @test result["termination_status"] == OPTIMAL
        @test isapprox(result["solution"]["node"]["1"]["h"], 10.0, rtol=1.0e-3)
        @test isapprox(result["solution"]["pipe"]["1"]["q"], 0.9, rtol=1.0e-3)
    end

    @testset "Head Loss (PRV)" begin
        network = WaterModels.parse_file("../test/data/epanet/snapshot/prv-hw-lps.inp")

        wm = instantiate_model(network, NLPWaterModel, build_wf)
        result = WaterModels.optimize_model!(wm, optimizer=_make_juniper(wm, ipopt))
        @test result["termination_status"] == LOCALLY_SOLVED
        @test isapprox(result["solution"]["node"]["1"]["h"], 10.0, rtol=1.0e-3)
        @test isapprox(result["solution"]["node"]["2"]["p"], 2.38, rtol=1.0e-3)
        @test isapprox(result["solution"]["node"]["3"]["h"], 2.00, rtol=1.0e-3)
        @test isapprox(result["solution"]["node"]["3"]["p"], 1.00, rtol=1.0e-3)
        @test isapprox(result["solution"]["pipe"]["1"]["q"], 2.0, rtol=1.0e-3)

        wm = instantiate_model(network, MICPRWaterModel, build_wf)
        result = WaterModels.optimize_model!(wm, optimizer=_make_juniper(wm, ipopt))
        @test result["termination_status"] == LOCALLY_SOLVED
        @test isapprox(result["solution"]["node"]["1"]["h"], 10.0, rtol=1.0e-3)
        @test isapprox(result["solution"]["node"]["3"]["h"], 2.00, rtol=1.0e-3)
        @test isapprox(result["solution"]["node"]["3"]["p"], 1.00, rtol=1.0e-3)
        @test isapprox(result["solution"]["pipe"]["1"]["q"], 2.0, rtol=1.0e-3)

        result = run_wf(network, MILPWaterModel, cbc, ext=Dict(:pipe_breakpoints=>3))
        @test result["termination_status"] == OPTIMAL
        @test isapprox(result["solution"]["node"]["1"]["h"], 10.0, rtol=1.0e-3)
        @test isapprox(result["solution"]["node"]["2"]["p"], 2.38, rtol=1.0e-3)
        @test isapprox(result["solution"]["node"]["3"]["h"], 2.00, rtol=1.0e-3)
        @test isapprox(result["solution"]["node"]["3"]["p"], 1.00, rtol=1.0e-3)
        @test isapprox(result["solution"]["pipe"]["1"]["q"], 2.0, rtol=1.0e-3)

        result = run_wf(network, MILPRWaterModel, cbc, ext=Dict(:pipe_breakpoints=>3))
        @test result["termination_status"] == OPTIMAL
        @test isapprox(result["solution"]["node"]["1"]["h"], 10.0, rtol=1.0e-3)
        @test isapprox(result["solution"]["node"]["3"]["h"], 2.00, rtol=1.0e-3)
        @test isapprox(result["solution"]["node"]["3"]["p"], 1.00, rtol=1.0e-3)
        @test isapprox(result["solution"]["pipe"]["1"]["q"], 2.0, rtol=1.0e-3)
    end

    @testset "Head Gain (Pump)" begin
        network = WaterModels.parse_file("../test/data/epanet/snapshot/pump-hw-lps.inp")

        wm = instantiate_model(network, NLPWaterModel, build_wf)
        result = WaterModels.optimize_model!(wm, optimizer=_make_juniper(wm, ipopt))
        @test result["termination_status"] == LOCALLY_SOLVED
        @test isapprox(result["solution"]["node"]["1"]["h"], 10.0, rtol=1.0e-3)
        @test isapprox(result["solution"]["node"]["2"]["h"], 98.98, rtol=1.0e-3)
        @test isapprox(result["solution"]["pump"]["1"]["status"], 1.0, atol=1.0e-3)

        wm = instantiate_model(network, MICPRWaterModel, build_wf)
        result = WaterModels.optimize_model!(wm, optimizer=_make_juniper(wm, ipopt))
        @test result["termination_status"] == LOCALLY_SOLVED
        @test isapprox(result["solution"]["node"]["1"]["h"], 10.0, rtol=1.0e-3)
        @test isapprox(result["solution"]["pump"]["1"]["status"], 1.0, atol=1.0e-3)

        result = run_wf(network, MILPWaterModel, cbc, ext=Dict(:pump_breakpoints=>4))
        @test result["termination_status"] == OPTIMAL
        @test isapprox(result["solution"]["node"]["1"]["h"], 10.0, rtol=1.0e-3)
        @test isapprox(result["solution"]["node"]["2"]["h"], 98.98, rtol=1.0e-1)
        @test isapprox(result["solution"]["pump"]["1"]["status"], 1.0, atol=1.0e-3)

        result = run_wf(network, MILPRWaterModel, cbc, ext=Dict(:pump_breakpoints=>3))
        @test result["termination_status"] == OPTIMAL
        @test isapprox(result["solution"]["node"]["1"]["h"], 10.0, rtol=1.0e-3)
        @test isapprox(result["solution"]["pump"]["1"]["status"], 1.0, atol=1.0e-3)
    end

    @testset "Hazen-Williams Head Loss (Tank)" begin
        network = parse_file("../test/data/epanet/snapshot/tank-hw-lps.inp")

        wm = instantiate_model(network, NLPWaterModel, build_wf)
        result = WaterModels.optimize_model!(wm, optimizer=_make_juniper(wm, ipopt))
        @test result["termination_status"] == LOCALLY_SOLVED
        @test isapprox(result["solution"]["node"]["1"]["h"], 20.0, rtol=1.0e-3)
        @test isapprox(result["solution"]["node"]["2"]["h"], 17.89, rtol=1.0e-3)
        @test isapprox(result["solution"]["pipe"]["1"]["q"], 1.0, rtol=1.0e-3)

        wm = instantiate_model(network, MICPRWaterModel, build_wf)
        result = WaterModels.optimize_model!(wm, optimizer=_make_juniper(wm, ipopt))
        @test result["termination_status"] == LOCALLY_SOLVED
        @test isapprox(result["solution"]["node"]["1"]["h"], 20.0, rtol=1.0e-3)
        @test result["solution"]["node"]["2"]["h"] <= 17.89
        @test isapprox(result["solution"]["pipe"]["1"]["q"], 1.0, rtol=1.0e-3)

        result = run_wf(network, MILPWaterModel, cbc, ext=Dict(:pipe_breakpoints=>4))
        @test result["termination_status"] == OPTIMAL
        @test isapprox(result["solution"]["node"]["1"]["h"], 20.0, rtol=1.0e-3)
        @test isapprox(result["solution"]["node"]["2"]["h"], 17.89, rtol=1.0e-1)
        @test isapprox(result["solution"]["pipe"]["1"]["q"], 1.0, rtol=1.0e-3)

        result = run_wf(network, MILPRWaterModel, cbc, ext=Dict(:pipe_breakpoints=>3))
        @test result["termination_status"] == OPTIMAL
        @test isapprox(result["solution"]["node"]["1"]["h"], 20.0, rtol=1.0e-3)
        @test isapprox(result["solution"]["pipe"]["1"]["q"], 1.0, rtol=1.0e-3)
    end
end

@testset "Water Flow Problems (Multinetwork)" begin
    # Pipe, check valve, and shutoff valve tests will have the same physical solutions.
    for (component, name) in Dict("Pipe"=>"pipe", "Check Valve"=>"check_valve", "Shutoff Valve"=>"shutoff_valve")
        @testset "Hazen-Williams Head Loss ($(component))" begin
            network = WaterModels.parse_file("../test/data/epanet/multinetwork/$(name)-hw-lps.inp")
            network = WaterModels.make_multinetwork(network)

            wm = instantiate_model(network, NLPWaterModel, build_mn_wf)
            result = WaterModels.optimize_model!(wm, optimizer=_make_juniper(wm, ipopt))
            @test result["termination_status"] == LOCALLY_SOLVED
            @test isapprox(result["solution"]["nw"]["2"]["node"]["1"]["h"], 10.0, rtol=1.0e-3)
            @test isapprox(result["solution"]["nw"]["1"]["node"]["2"]["p"], 7.89, rtol=1.0e-3)
            @test isapprox(result["solution"]["nw"]["2"]["node"]["2"]["h"], 9.42, rtol=1.0e-3)
            @test isapprox(result["solution"]["nw"]["3"]["node"]["2"]["h"], 9.84, rtol=1.0e-3)
            @test isapprox(result["solution"]["nw"]["1"][name]["1"]["q"], 1.0, rtol=1.0e-3)
            @test isapprox(result["solution"]["nw"]["2"][name]["1"]["q"], 0.5, rtol=1.0e-3)
            @test isapprox(result["solution"]["nw"]["3"][name]["1"]["q"], 0.25, rtol=1.0e-3)

            wm = instantiate_model(network, MICPRWaterModel, build_mn_wf)
            result = WaterModels.optimize_model!(wm, optimizer=_make_juniper(wm, ipopt))
            @test result["termination_status"] == LOCALLY_SOLVED
            @test isapprox(result["solution"]["nw"]["2"]["node"]["1"]["h"], 10.0, rtol=1.0e-3)
            @test isapprox(result["solution"]["nw"]["1"][name]["1"]["q"], 1.0, rtol=1.0e-3)
            @test isapprox(result["solution"]["nw"]["2"][name]["1"]["q"], 0.5, rtol=1.0e-3)
            @test isapprox(result["solution"]["nw"]["3"][name]["1"]["q"], 0.25, rtol=1.0e-3)

            result = run_mn_wf(network, MILPWaterModel, cbc, ext=Dict(:pipe_breakpoints=>3))
            @test result["termination_status"] == OPTIMAL
            @test isapprox(result["solution"]["nw"]["2"]["node"]["1"]["h"], 10.0, rtol=1.0e-3)
            @test isapprox(result["solution"]["nw"]["1"]["node"]["2"]["p"], 7.89, rtol=1.0e-3)
            @test isapprox(result["solution"]["nw"]["2"]["node"]["2"]["h"], 9.42, rtol=1.0e-3)
            @test isapprox(result["solution"]["nw"]["3"]["node"]["2"]["h"], 9.84, rtol=1.0e-3)
            @test isapprox(result["solution"]["nw"]["1"][name]["1"]["q"], 1.0, rtol=1.0e-3)
            @test isapprox(result["solution"]["nw"]["2"][name]["1"]["q"], 0.5, rtol=1.0e-3)
            @test isapprox(result["solution"]["nw"]["3"][name]["1"]["q"], 0.25, rtol=1.0e-3)

            result = run_mn_wf(network, MILPRWaterModel, cbc, ext=Dict(:pipe_breakpoints=>3))
            @test result["termination_status"] == OPTIMAL
            @test isapprox(result["solution"]["nw"]["2"]["node"]["1"]["h"], 10.0, rtol=1.0e-3)
            @test isapprox(result["solution"]["nw"]["1"][name]["1"]["q"], 1.0, rtol=1.0e-3)
            @test isapprox(result["solution"]["nw"]["2"][name]["1"]["q"], 0.5, rtol=1.0e-3)
            @test isapprox(result["solution"]["nw"]["3"][name]["1"]["q"], 0.25, rtol=1.0e-3)
        end
    end

    @testset "Hazen-Williams Head Loss (Negative Demand)" begin
        network = WaterModels.parse_file("../test/data/epanet/multinetwork/negative_demand-hw-lps.inp")
        network = WaterModels.make_multinetwork(network)

        wm = instantiate_model(network, NLPWaterModel, build_mn_wf)
        result = WaterModels.optimize_model!(wm, optimizer=_make_juniper(wm, ipopt))
        @test result["termination_status"] == LOCALLY_SOLVED
        @test isapprox(result["solution"]["nw"]["1"]["node"]["1"]["h"], 10.0, rtol=1.0e-3)
        @test isapprox(result["solution"]["nw"]["2"]["node"]["2"]["p"], 8.27, rtol=1.0e-3)
        @test isapprox(result["solution"]["nw"]["3"]["node"]["3"]["p"], 8.29, rtol=1.0e-3)
        @test isapprox(result["solution"]["nw"]["1"]["pipe"]["1"]["q"], 0.9, rtol=1.0e-3)

        wm = instantiate_model(network, MICPRWaterModel, build_mn_wf)
        result = WaterModels.optimize_model!(wm, optimizer=_make_juniper(wm, ipopt))
        @test result["termination_status"] == LOCALLY_SOLVED
        @test isapprox(result["solution"]["nw"]["1"]["node"]["1"]["h"], 10.0, rtol=1.0e-3)
        @test isapprox(result["solution"]["nw"]["3"]["pipe"]["1"]["q"], 0.9, rtol=1.0e-3)

        result = run_mn_wf(network, MILPWaterModel, cbc, ext=Dict(:pipe_breakpoints=>4))
        @test result["termination_status"] == OPTIMAL
        @test isapprox(result["solution"]["nw"]["1"]["node"]["1"]["h"], 10.0, rtol=1.0e-3)
        @test isapprox(result["solution"]["nw"]["2"]["node"]["2"]["p"], 8.27, rtol=1.0e-2)
        @test isapprox(result["solution"]["nw"]["3"]["node"]["3"]["p"], 8.29, rtol=1.0e-3)
        @test isapprox(result["solution"]["nw"]["1"]["pipe"]["1"]["q"], 0.9, rtol=1.0e-3)

        result = run_mn_wf(network, MILPRWaterModel, cbc, ext=Dict(:pipe_breakpoints=>3))
        @test result["termination_status"] == OPTIMAL
        @test isapprox(result["solution"]["nw"]["1"]["node"]["1"]["h"], 10.0, rtol=1.0e-3)
        @test isapprox(result["solution"]["nw"]["2"]["pipe"]["1"]["q"], 0.9, rtol=1.0e-3)
    end

    @testset "Head Loss (PRV)" begin
        network = WaterModels.parse_file("../test/data/epanet/multinetwork/prv-hw-lps.inp")
        network = WaterModels.make_multinetwork(network)

        wm = instantiate_model(network, NLPWaterModel, build_mn_wf)
        result = WaterModels.optimize_model!(wm, optimizer=_make_juniper(wm, ipopt))
        @test result["termination_status"] == LOCALLY_SOLVED
        @test isapprox(result["solution"]["nw"]["1"]["node"]["2"]["h"], 2.38, rtol=1.0e-3)
        @test isapprox(result["solution"]["nw"]["2"]["node"]["2"]["h"], 7.89, rtol=1.0e-3)
        @test isapprox(result["solution"]["nw"]["3"]["node"]["2"]["h"], 9.42, rtol=1.0e-3)
        @test isapprox(result["solution"]["nw"]["2"]["node"]["3"]["h"], 2.00, rtol=1.0e-3)
        @test isapprox(result["solution"]["nw"]["2"]["node"]["1"]["h"], 10.0, rtol=1.0e-3)
        @test isapprox(result["solution"]["nw"]["1"]["pipe"]["1"]["q"], 2.00, rtol=1.0e-3)
        @test isapprox(result["solution"]["nw"]["3"]["pipe"]["1"]["q"], 0.50, rtol=1.0e-3)
        @test isapprox(result["solution"]["nw"]["1"]["pressure_reducing_valve"]["2"]["q"], 1.00, rtol=1.0e-3)
        @test isapprox(result["solution"]["nw"]["3"]["pressure_reducing_valve"]["2"]["q"], 0.25, rtol=1.0e-3)

        wm = instantiate_model(network, MICPRWaterModel, build_mn_wf)
        result = WaterModels.optimize_model!(wm, optimizer=_make_juniper(wm, ipopt))
        @test result["termination_status"] == LOCALLY_SOLVED
        @test result["solution"]["nw"]["1"]["node"]["2"]["h"] <= 2.39
        @test result["solution"]["nw"]["2"]["node"]["2"]["h"] <= 7.90
        @test result["solution"]["nw"]["3"]["node"]["2"]["h"] <= 9.42
        @test isapprox(result["solution"]["nw"]["2"]["node"]["3"]["h"], 2.00, rtol=1.0e-3)
        @test isapprox(result["solution"]["nw"]["2"]["node"]["1"]["h"], 10.0, rtol=1.0e-3)
        @test isapprox(result["solution"]["nw"]["1"]["pipe"]["1"]["q"], 2.00, rtol=1.0e-3)
        @test isapprox(result["solution"]["nw"]["3"]["pipe"]["1"]["q"], 0.50, rtol=1.0e-3)
        @test isapprox(result["solution"]["nw"]["1"]["pressure_reducing_valve"]["2"]["q"], 1.00, rtol=1.0e-3)
        @test isapprox(result["solution"]["nw"]["3"]["pressure_reducing_valve"]["2"]["q"], 0.25, rtol=1.0e-3)

        result = run_mn_wf(network, MILPWaterModel, cbc, ext=Dict(:pipe_breakpoints=>3))
        @test result["termination_status"] == OPTIMAL
        @test isapprox(result["solution"]["nw"]["1"]["node"]["2"]["h"], 2.38, rtol=1.0e-3)
        @test isapprox(result["solution"]["nw"]["2"]["node"]["2"]["h"], 7.89, rtol=1.0e-3)
        @test isapprox(result["solution"]["nw"]["3"]["node"]["2"]["h"], 9.42, rtol=1.0e-3)
        @test isapprox(result["solution"]["nw"]["2"]["node"]["3"]["h"], 2.00, rtol=1.0e-3)
        @test isapprox(result["solution"]["nw"]["2"]["node"]["1"]["h"], 10.0, rtol=1.0e-3)
        @test isapprox(result["solution"]["nw"]["1"]["pipe"]["1"]["q"], 2.00, rtol=1.0e-3)
        @test isapprox(result["solution"]["nw"]["3"]["pipe"]["1"]["q"], 0.50, rtol=1.0e-3)
        @test isapprox(result["solution"]["nw"]["1"]["pressure_reducing_valve"]["2"]["q"], 1.00, rtol=1.0e-3)
        @test isapprox(result["solution"]["nw"]["3"]["pressure_reducing_valve"]["2"]["q"], 0.25, rtol=1.0e-3)

        result = run_mn_wf(network, MILPRWaterModel, cbc, ext=Dict(:pipe_breakpoints=>3))
        @test result["termination_status"] == OPTIMAL
        @test result["solution"]["nw"]["1"]["node"]["2"]["h"] <= 2.39
        @test result["solution"]["nw"]["2"]["node"]["2"]["h"] <= 7.90
        @test result["solution"]["nw"]["3"]["node"]["2"]["h"] <= 9.42
        @test isapprox(result["solution"]["nw"]["2"]["node"]["3"]["h"], 2.00, rtol=1.0e-3)
        @test isapprox(result["solution"]["nw"]["2"]["node"]["1"]["h"], 10.0, rtol=1.0e-3)
        @test isapprox(result["solution"]["nw"]["1"]["pipe"]["1"]["q"], 2.00, rtol=1.0e-3)
        @test isapprox(result["solution"]["nw"]["3"]["pipe"]["1"]["q"], 0.50, rtol=1.0e-3)
        @test isapprox(result["solution"]["nw"]["1"]["pressure_reducing_valve"]["2"]["q"], 1.00, rtol=1.0e-3)
        @test isapprox(result["solution"]["nw"]["3"]["pressure_reducing_valve"]["2"]["q"], 0.25, rtol=1.0e-3)
    end

    @testset "Head Gain (Pump)" begin
        network = WaterModels.parse_file("../test/data/epanet/multinetwork/pump-hw-lps.inp")
        network = WaterModels.make_multinetwork(network)

        wm = instantiate_model(network, NLPWaterModel, build_mn_wf)
        result = WaterModels.optimize_model!(wm, optimizer=_make_juniper(wm, ipopt))
        @test result["termination_status"] == LOCALLY_SOLVED
        @test isapprox(result["solution"]["nw"]["1"]["pump"]["1"]["q"], 0.125, rtol=1.0e-3)
        @test isapprox(result["solution"]["nw"]["3"]["pump"]["1"]["q"], 0.03125, rtol=1.0e-3)
        @test isapprox(result["solution"]["nw"]["1"]["pump"]["1"]["status"], 1.0, atol=1.0e-3)
        @test isapprox(result["solution"]["nw"]["1"]["pump"]["1"]["g"], 88.98, rtol=1.0e-3)
        @test isapprox(result["solution"]["nw"]["3"]["pump"]["1"]["g"], 99.59, rtol=1.0e-2)

        wm = instantiate_model(network, MICPRWaterModel, build_mn_wf)
        result = WaterModels.optimize_model!(wm, optimizer=_make_juniper(wm, ipopt))
        @test result["termination_status"] == LOCALLY_SOLVED
        @test isapprox(result["solution"]["nw"]["1"]["pump"]["1"]["q"], 0.125, rtol=1.0e-3)
        @test isapprox(result["solution"]["nw"]["3"]["pump"]["1"]["q"], 0.03125, rtol=1.0e-3)
        @test isapprox(result["solution"]["nw"]["1"]["pump"]["1"]["status"], 1.0, atol=1.0e-3)
        @test result["solution"]["nw"]["1"]["pump"]["1"]["g"] <= 88.99
        @test result["solution"]["nw"]["3"]["pump"]["1"]["g"] <= 99.60

        result = run_mn_wf(network, MILPWaterModel, cbc, ext=Dict(:pump_breakpoints=>5))
        @test result["termination_status"] == OPTIMAL
        @test isapprox(result["solution"]["nw"]["1"]["pump"]["1"]["q"], 0.125, rtol=1.0e-3)
        @test isapprox(result["solution"]["nw"]["3"]["pump"]["1"]["q"], 0.03125, rtol=1.0e-3)
        @test isapprox(result["solution"]["nw"]["1"]["pump"]["1"]["status"], 1.0, atol=1.0e-3)
        @test isapprox(result["solution"]["nw"]["1"]["pump"]["1"]["g"], 88.98, rtol=1.0e-1)
        @test isapprox(result["solution"]["nw"]["3"]["pump"]["1"]["g"], 99.59, rtol=1.0e-1)

        result = run_mn_wf(network, MILPRWaterModel, cbc, ext=Dict(:pump_breakpoints=>3))
        @test result["termination_status"] == OPTIMAL
        @test isapprox(result["solution"]["nw"]["1"]["pump"]["1"]["q"], 0.125, rtol=1.0e-3)
        @test isapprox(result["solution"]["nw"]["3"]["pump"]["1"]["q"], 0.03125, rtol=1.0e-3)
        @test isapprox(result["solution"]["nw"]["1"]["pump"]["1"]["status"], 1.0, atol=1.0e-3)
        @test result["solution"]["nw"]["1"]["pump"]["1"]["g"] <= 88.99
        @test result["solution"]["nw"]["3"]["pump"]["1"]["g"] <= 99.60
    end

    @testset "Hazen-Williams Head Loss (Tank)" begin
        network = WaterModels.parse_file("../test/data/epanet/multinetwork/tank-hw-lps.inp")
        network = WaterModels.make_multinetwork(network)

        wm = instantiate_model(network, NLPWaterModel, build_mn_wf)
        result = WaterModels.optimize_model!(wm, optimizer=_make_juniper(wm, ipopt))
        @test result["termination_status"] == LOCALLY_SOLVED
        @test isapprox(result["solution"]["nw"]["1"]["node"]["1"]["h"], 60.00, rtol=1.0e-3)
        @test isapprox(result["solution"]["nw"]["1"]["node"]["2"]["h"], 59.42, rtol=1.0e-3)
        @test isapprox(result["solution"]["nw"]["3"]["node"]["1"]["h"], 25.62, rtol=1.0e-3)
        @test isapprox(result["solution"]["nw"]["3"]["node"]["2"]["h"], 25.58, rtol=1.0e-3)
        @test isapprox(result["solution"]["nw"]["1"]["pipe"]["1"]["q"], 0.500, rtol=1.0e-3)
        @test isapprox(result["solution"]["nw"]["3"]["pipe"]["1"]["q"], 0.125, rtol=1.0e-3)

        wm = instantiate_model(network, MICPRWaterModel, build_mn_wf)
        result = WaterModels.optimize_model!(wm, optimizer=_make_juniper(wm, ipopt))
        @test result["termination_status"] == LOCALLY_SOLVED
        @test isapprox(result["solution"]["nw"]["1"]["node"]["1"]["h"], 60.00, rtol=1.0e-3)
        @test result["solution"]["nw"]["1"]["node"]["2"]["h"] <= 59.42
        @test isapprox(result["solution"]["nw"]["3"]["node"]["1"]["h"], 25.62, rtol=1.0e-3)
        @test result["solution"]["nw"]["3"]["node"]["2"]["h"] <= 25.58
        @test isapprox(result["solution"]["nw"]["1"]["pipe"]["1"]["q"], 0.500, rtol=1.0e-3)
        @test isapprox(result["solution"]["nw"]["3"]["pipe"]["1"]["q"], 0.125, rtol=1.0e-3)

        result = run_mn_wf(network, MILPWaterModel, cbc, ext=Dict(:pipe_breakpoints=>4))
        @test result["termination_status"] == OPTIMAL
        @test isapprox(result["solution"]["nw"]["1"]["node"]["1"]["h"], 60.00, rtol=1.0e-3)
        @test isapprox(result["solution"]["nw"]["1"]["node"]["2"]["h"], 59.42, rtol=1.0e-2)
        @test isapprox(result["solution"]["nw"]["3"]["node"]["1"]["h"], 25.62, rtol=1.0e-3)
        @test isapprox(result["solution"]["nw"]["3"]["node"]["2"]["h"], 25.58, rtol=1.0e-2)
        @test isapprox(result["solution"]["nw"]["1"]["pipe"]["1"]["q"], 0.500, rtol=1.0e-3)
        @test isapprox(result["solution"]["nw"]["3"]["pipe"]["1"]["q"], 0.125, rtol=1.0e-3)

        result = run_mn_wf(network, MILPRWaterModel, cbc, ext=Dict(:pipe_breakpoints=>3))
        @test result["termination_status"] == OPTIMAL
        @test isapprox(result["solution"]["nw"]["1"]["node"]["1"]["h"], 60.00, rtol=1.0e-3)
        @test isapprox(result["solution"]["nw"]["3"]["node"]["1"]["h"], 25.62, rtol=1.0e-3)
        @test isapprox(result["solution"]["nw"]["1"]["pipe"]["1"]["q"], 0.500, rtol=1.0e-3)
        @test isapprox(result["solution"]["nw"]["3"]["pipe"]["1"]["q"], 0.125, rtol=1.0e-3)
    end
end
