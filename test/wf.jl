@testset "Water Flow Problems" begin
    @testset "2PRVs network, MICP formulation." begin
        data = WaterModels.parse_file("../test/data/epanet/2PRVs.inp")
        wm = instantiate_model(data, MICPWaterModel, WaterModels.build_wf)
        f = Juniper.register(head_loss_args(wm)..., autodiff=false)
        juniper = JuMP.optimizer_with_attributes(Juniper.Optimizer,
            "nl_solver"=>ipopt, "registered_functions"=>[f], "log_levels"=>[])
        solution = WaterModels.optimize_model!(wm, optimizer=juniper)
        @test solution["termination_status"] == LOCALLY_SOLVED
    end

    @testset "Balerma network, MICP formulation." begin
        data = WaterModels.parse_file("../test/data/epanet/balerma.inp")
        modifications = WaterModels.parse_file("../test/data/json/balerma.json")
        InfrastructureModels.update_data!(data, modifications)
        wm = instantiate_model(data, MICPWaterModel, WaterModels.build_wf)
        f = Juniper.register(head_loss_args(wm)..., autodiff=false)
        juniper = JuMP.optimizer_with_attributes(Juniper.Optimizer, "nl_solver"=>ipopt,
            "registered_functions"=>[f], "log_levels"=>[])

        # TODO: The below takes too long to execute.
        #solution = WaterModels.optimize_model!(wm, optimizer=juniper)
        #@test solution["termination_status"] == LOCALLY_SOLVED
    end

    @testset "Example 1 network (with tanks and pumps), MICP formulation." begin
        data = WaterModels.parse_file("../test/data/epanet/example_1-sp.inp")
        wm = instantiate_model(data, MICPWaterModel, WaterModels.build_wf)
        f = Juniper.register(head_loss_args(wm)..., autodiff=false)
        juniper = JuMP.optimizer_with_attributes(Juniper.Optimizer,
            "nl_solver"=>ipopt, "registered_functions"=>[f], "log_levels"=>[])
        solution = WaterModels.optimize_model!(wm, optimizer=juniper)
        @test solution["termination_status"] == LOCALLY_SOLVED
    end

    @testset "Example 1 network (with tanks and pumps), multinetwork MICP formulation." begin
        data = WaterModels.parse_file("../test/data/epanet/example_1.inp")
        mn_data = WaterModels.make_multinetwork(data)
        wm = instantiate_model(mn_data, MICPWaterModel, WaterModels.build_mn_wf)
        f = Juniper.register(head_loss_args(wm)..., autodiff=false)
        juniper = JuMP.optimizer_with_attributes(Juniper.Optimizer,
            "nl_solver"=>ipopt, "registered_functions"=>[f], "log_levels"=>[])
        solution = WaterModels.optimize_model!(wm, optimizer=juniper)
        @test solution["termination_status"] == LOCALLY_SOLVED
    end

    @testset "Richmond (single time) network, MICP formulation." begin
        data = WaterModels.parse_file("../test/data/epanet/richmond-skeleton-sp.inp")
        wm = instantiate_model(data, MICPWaterModel, WaterModels.build_wf)
        f = Juniper.register(head_loss_args(wm)..., autodiff=false)
        juniper = JuMP.optimizer_with_attributes(Juniper.Optimizer,
            "nl_solver"=>ipopt, "registered_functions"=>[f], "log_levels"=>[])
        # TODO: The below takes too long to execute.
        #solution = WaterModels.optimize_model!(wm, optimizer=juniper)
        #@test solution["termination_status"] == LOCALLY_SOLVED
    end

    @testset "Richmond network, multinetwork MICP formulation." begin
        data = WaterModels.parse_file("../test/data/epanet/richmond-skeleton.inp")
        mn_data = WaterModels.make_multinetwork(data)
        wm = instantiate_model(mn_data, MICPWaterModel, WaterModels.build_mn_wf)
        f = Juniper.register(head_loss_args(wm)..., autodiff=false)
        juniper = JuMP.optimizer_with_attributes(Juniper.Optimizer,
            "nl_solver"=>ipopt, "registered_functions"=>[f], "log_levels"=>[])
        # TODO: The below takes too long to execute.
        #solution = WaterModels.optimize_model!(wm, optimizer=juniper)
        #@test solution["termination_status"] == LOCALLY_SOLVED
    end

    @testset "Shamir network, MICP formulation." begin
        data = WaterModels.parse_file("../test/data/epanet/shamir.inp")
        wm = instantiate_model(data, MICPWaterModel, WaterModels.build_wf)
        f = Juniper.register(head_loss_args(wm)..., autodiff=false)
        juniper = JuMP.optimizer_with_attributes(Juniper.Optimizer,
            "nl_solver"=>ipopt, "registered_functions"=>[f], "log_levels"=>[])
        solution = WaterModels.optimize_model!(wm, optimizer=juniper)
        @test solution["termination_status"] == LOCALLY_SOLVED
    end

    @testset "Shamir network, multinetwork MICP formulation." begin
        data = WaterModels.parse_file("../test/data/epanet/shamir-ts.inp")
        mn_data = WaterModels.make_multinetwork(data)
        wm = instantiate_model(mn_data, MICPWaterModel, WaterModels.build_mn_wf)
        f = Juniper.register(head_loss_args(wm)..., autodiff=false)
        juniper = JuMP.optimizer_with_attributes(Juniper.Optimizer,
            "nl_solver"=>ipopt, "registered_functions"=>[f], "log_levels"=>[])
        solution = WaterModels.optimize_model!(wm, optimizer=juniper)
        @test solution["termination_status"] == LOCALLY_SOLVED
    end

    @testset "Shamir network (with tank), multinetwork MICP formulation." begin
        data = WaterModels.parse_file("../test/data/epanet/shamir-ts-tank.inp")
        mn_data = WaterModels.make_multinetwork(data)
        wm = instantiate_model(mn_data, MICPWaterModel, WaterModels.build_mn_wf)
        f = Juniper.register(head_loss_args(wm)..., autodiff=false)
        juniper = JuMP.optimizer_with_attributes(Juniper.Optimizer,
            "nl_solver"=>ipopt, "registered_functions"=>[f], "log_levels"=>[])
        solution = WaterModels.optimize_model!(wm, optimizer=juniper)
        @test solution["termination_status"] == LOCALLY_SOLVED
    end

    @testset "2PRVs network, MILP formulation." begin
        ext = Dict(:pipe_breakpoints=>25)
        solution = solve_wf("../test/data/epanet/2PRVs.inp", MILPWaterModel, cbc, ext=ext)

        @test solution["termination_status"] == OPTIMAL
        @test isapprox(solution["solution"]["link"]["1"]["q"], 0.012618, rtol=1.0e-3)
        @test isapprox(solution["solution"]["link"]["4"]["q"], 0.0, atol=1.0e-4)
        @test isapprox(solution["solution"]["link"]["5"]["q"], 0.012618, rtol=1.0e-3)
        @test isapprox(solution["solution"]["node"]["1"]["h"], 45.720001, rtol=1.0e-3)
        @test isapprox(solution["solution"]["node"]["3"]["h"], 35.166771, rtol=1.0e-3)
    end

    @testset "Balerma network, MILP formulation." begin
        ext = Dict(:pipe_breakpoints=>10, :pump_breakpoints=>10)
        data = WaterModels.parse_file("../test/data/epanet/balerma.inp")
        modifications = WaterModels.parse_file("../test/data/json/balerma.json")
        InfrastructureModels.update_data!(data, modifications)
        wm = instantiate_model(data, MILPWaterModel, WaterModels.build_wf, ext=ext)

        ## TODO: The below takes too long to solve with cbc.
        #solution = WaterModels.optimize_model!(wm, optimizer=cbc)
        #@test solution["termination_status"] == OPTIMAL
        #@test isapprox(solution["solution"]["link"]["1"]["q"], -0.002498, rtol=1.0e-3)
        #@test isapprox(solution["solution"]["link"]["429"]["q"], 0.000864, rtol=1.0e-3)
        #@test isapprox(solution["solution"]["link"]["540"]["q"], -0.008713, rtol=1.0e-3)
        #@test isapprox(solution["solution"]["node"]["191"]["h"], 118.006958, rtol=1.0e-3)
        #@test isapprox(solution["solution"]["node"]["381"]["h"], 86.70694, rtol=1.0e-3)
        #@test isapprox(solution["solution"]["node"]["127001"]["h"], 84.481102, rtol=1.0e-3)
    end

    @testset "Example 1 network (with tanks and pump), MILP formulation." begin
        data = WaterModels.parse_file("../test/data/epanet/example_1-sp.inp")
        ext = Dict(:pipe_breakpoints=>25, :pump_breakpoints=>25)
        wm = instantiate_model(data, MILPWaterModel, WaterModels.build_wf, ext=ext)

        # TODO: Uncomment once the CBC bug is fixed.
        #solution = _IM.optimize_model!(wm, optimizer=cbc)
        #@test solution["termination_status"] == OPTIMAL
        #@test isapprox(solution["solution"]["link"]["9"]["q"], 0.117737, rtol=1.0e-2)
        #@test isapprox(solution["solution"]["link"]["110"]["q"], -0.048338, rtol=1.0e-2)
        #@test isapprox(solution["solution"]["link"]["122"]["q"], 0.003734, rtol=1.0e-2)
        #@test isapprox(solution["solution"]["node"]["9"]["h"], 243.839996, rtol=1.0e-2)
        #@test isapprox(solution["solution"]["node"]["10"]["h"], 306.125092, rtol=1.0e-2)
        #@test isapprox(solution["solution"]["node"]["23"]["h"], 295.243073, rtol=1.0e-2)
    end

    @testset "Example 1 network (with tanks and links), multinetwork MILP formulation." begin
        ext = Dict(:pipe_breakpoints=>12, :pump_breakpoints=>12)
        data = WaterModels.parse_file("../test/data/epanet/example_1.inp")
        mn_data = WaterModels.make_multinetwork(data)
        wm = instantiate_model(mn_data, MILPWaterModel, WaterModels.build_mn_wf, ext=ext)
        solution = _IM.optimize_model!(wm, optimizer=cbc)
        @test solution["termination_status"] == OPTIMAL
    end

    @testset "Richmond (single time) network, MILP formulation." begin
        ext = Dict(:pipe_breakpoints=>20, :pump_breakpoints=>20)
        data = WaterModels.parse_file("../test/data/epanet/richmond-skeleton-sp.inp")
        wm = instantiate_model(data, MILPWaterModel, WaterModels.build_wf, ext=ext)
        solution = _IM.optimize_model!(wm, optimizer=cbc)
        @test solution["termination_status"] == OPTIMAL
    end

    @testset "Shamir network, MILP formulation." begin
        ext = Dict(:pipe_breakpoints=>10, :pump_breakpoints=>10)
        solution = solve_wf("../test/data/epanet/shamir.inp", MILPWaterModel, cbc, ext=ext)

        @test solution["termination_status"] == OPTIMAL
        @test isapprox(solution["solution"]["link"]["2"]["q"], 0.093565, rtol=1.0e-2)
        @test isapprox(solution["solution"]["link"]["6"]["q"], 0.055710, rtol=1.0e-2)
        @test isapprox(solution["solution"]["node"]["2"]["h"], 203.247650, rtol=1.0e-2)
        @test isapprox(solution["solution"]["node"]["6"]["h"], 195.445953, rtol=1.0e-2)
        @test isapprox(solution["solution"]["reservoir"]["1"]["qr"], 0.31109, rtol=1.0e-2)
    end

    @testset "Shamir network, multinetwork MILP formulation." begin
        ext = Dict(:pipe_breakpoints=>10, :pump_breakpoints=>10)
        data = WaterModels.parse_file("../test/data/epanet/shamir-ts.inp")
        mn_data = WaterModels.make_multinetwork(data)
        wm = instantiate_model(mn_data, MILPWaterModel, WaterModels.build_mn_wf, ext=ext)
        solution = _IM.optimize_model!(wm, optimizer=cbc)

        @test solution["termination_status"] == OPTIMAL
        @test isapprox(solution["solution"]["nw"]["1"]["link"]["2"]["q"], 0.09356500, rtol=1.0e-2)
        @test isapprox(solution["solution"]["nw"]["1"]["link"]["6"]["q"], 0.05571000, rtol=1.0e-2)
        @test isapprox(solution["solution"]["nw"]["1"]["reservoir"]["1"]["qr"], 0.31109, rtol=1.0e-2)
        @test isapprox(solution["solution"]["nw"]["2"]["link"]["2"]["q"], 0.04678500, rtol=1.0e-2)
        @test isapprox(solution["solution"]["nw"]["2"]["link"]["6"]["q"], 0.02785300, rtol=1.0e-2)
        @test isapprox(solution["solution"]["nw"]["2"]["reservoir"]["1"]["qr"], 0.5*0.31109, rtol=1.0e-2)
        @test isapprox(solution["solution"]["nw"]["3"]["link"]["2"]["q"], 0.02339300, rtol=1.0e-2)
        @test isapprox(solution["solution"]["nw"]["3"]["link"]["6"]["q"], 0.01392600, rtol=1.0e-2)
        @test isapprox(solution["solution"]["nw"]["3"]["reservoir"]["1"]["qr"], 0.25*0.31109, rtol=1.0e-2)
    end

    @testset "Shamir network (with tank), multinetwork MILP formulation." begin
        ext = Dict(:pipe_breakpoints=>15, :pump_breakpoints=>15)
        data = WaterModels.parse_file("../test/data/epanet/shamir-ts-tank.inp")
        mn_data = WaterModels.make_multinetwork(data)
        wm = instantiate_model(mn_data, MILPWaterModel, WaterModels.build_mn_wf, ext=ext)
        solution = _IM.optimize_model!(wm, optimizer=cbc)

        @test solution["termination_status"] == OPTIMAL
        @test isapprox(solution["solution"]["nw"]["1"]["link"]["2"]["q"], 0.09356500, rtol=1.0e-2)
        @test isapprox(solution["solution"]["nw"]["1"]["link"]["6"]["q"], 0.05571000, rtol=1.0e-2)
        @test isapprox(solution["solution"]["nw"]["2"]["link"]["2"]["q"], 0.04678500, rtol=1.0e-2)
        @test isapprox(solution["solution"]["nw"]["2"]["link"]["6"]["q"], 0.02785300, rtol=1.0e-2)
        @test isapprox(solution["solution"]["nw"]["3"]["link"]["2"]["q"], 0.02339300, rtol=1.0e-2)
        @test isapprox(solution["solution"]["nw"]["3"]["link"]["6"]["q"], 0.01392600, rtol=1.0e-2)
        @test isapprox(solution["solution"]["nw"]["1"]["node"]["1"]["h"], 220.000000, rtol=1.0e-2)
        @test isapprox(solution["solution"]["nw"]["1"]["node"]["3"]["h"], 200.466599, rtol=1.0e-2)
        @test isapprox(solution["solution"]["nw"]["1"]["node"]["5"]["h"], 193.808380, rtol=1.0e-2)
        @test isapprox(solution["solution"]["nw"]["2"]["node"]["1"]["h"], 216.435196, rtol=1.0e-2)
        @test isapprox(solution["solution"]["nw"]["2"]["node"]["3"]["h"], 211.023941, rtol=1.0e-2)
        @test isapprox(solution["solution"]["nw"]["2"]["node"]["5"]["h"], 209.179306, rtol=1.0e-2)
        @test isapprox(solution["solution"]["nw"]["3"]["node"]["1"]["h"], 214.652771, rtol=1.0e-2)
        @test isapprox(solution["solution"]["nw"]["3"]["node"]["3"]["h"], 213.153824, rtol=1.0e-2)
        @test isapprox(solution["solution"]["nw"]["3"]["node"]["5"]["h"], 212.642853, rtol=1.0e-2)
    end

    @testset "2PRVs network, MILPR formulation." begin
        ext = Dict(:pipe_breakpoints=>5)
        solution = solve_wf("../test/data/epanet/2PRVs.inp", MILPRWaterModel, cbc, ext=ext)
        @test solution["termination_status"] == OPTIMAL
    end

    @testset "Balerma network, MILPR formulation." begin
        data = WaterModels.parse_file("../test/data/epanet/balerma.inp")
        modifications = WaterModels.parse_file("../test/data/json/balerma.json")
        InfrastructureModels.update_data!(data, modifications)
        ext = Dict(:pipe_breakpoints=>5, :pump_breakpoints=>5)
        wm = instantiate_model(data, MILPRWaterModel, WaterModels.build_wf, ext=ext)
        solution = WaterModels.optimize_model!(wm, optimizer=cbc)
        @test solution["termination_status"] == OPTIMAL
    end

    @testset "Example 1 network (with tanks and pumps), MILPR formulation." begin
        ext = Dict(:pipe_breakpoints=>5, :pump_breakpoints=>5)
        data = WaterModels.parse_file("../test/data/epanet/example_1-sp.inp")
        wm = instantiate_model(data, MILPRWaterModel, WaterModels.build_wf, ext=ext)
        juniper = JuMP.optimizer_with_attributes(Juniper.Optimizer, "nl_solver"=>ipopt, "log_levels"=>[])
        ## TODO: Uncomment this once CBC bug is fixed.
        #solution = WaterModels.optimize_model!(wm, optimizer=juniper)
        #@test solution["termination_status"] == LOCALLY_SOLVED
    end

    @testset "Example 1 network (with tanks and pumps), multinetwork MILPR formulation." begin
        ext = Dict(:pipe_breakpoints=>5, :pump_breakpoints=>5)
        data = WaterModels.parse_file("../test/data/epanet/example_1.inp")
        mn_data = WaterModels.make_multinetwork(data)
        wm = instantiate_model(mn_data, MILPRWaterModel, WaterModels.build_mn_wf, ext=ext)
        juniper = JuMP.optimizer_with_attributes(Juniper.Optimizer, "nl_solver"=>ipopt, "log_levels"=>[])
        solution = WaterModels.optimize_model!(wm, optimizer=juniper)
        @test solution["termination_status"] == LOCALLY_SOLVED
    end

    @testset "Richmond (single time) network, MILPR formulation." begin
        ext = Dict(:pipe_breakpoints=>5, :pump_breakpoints=>5)
        data = WaterModels.parse_file("../test/data/epanet/richmond-skeleton-sp.inp")
        wm = instantiate_model(data, MILPRWaterModel, WaterModels.build_wf, ext=ext)
        juniper = JuMP.optimizer_with_attributes(Juniper.Optimizer, "nl_solver"=>ipopt, "log_levels"=>[])
        # TODO: The below is stated to be infeasible.
        #solution = WaterModels.optimize_model!(wm, optimizer=juniper)
        #@test solution["termination_status"] == LOCALLY_SOLVED
    end

    @testset "Richmond network, multinetwork MILPR formulation." begin
        ext = Dict(:pipe_breakpoints=>5, :pump_breakpoints=>5)
        data = WaterModels.parse_file("../test/data/epanet/richmond-skeleton.inp")
        mn_data = WaterModels.make_multinetwork(data)
        wm = instantiate_model(mn_data, MILPRWaterModel, WaterModels.build_mn_wf, ext=ext)
        juniper = JuMP.optimizer_with_attributes(Juniper.Optimizer, "nl_solver"=>ipopt, "log_levels"=>[])
        # TODO: The below takes too long to execute.
        #solution = WaterModels.optimize_model!(wm, optimizer=juniper)
        #@test solution["termination_status"] == LOCALLY_SOLVED
    end

    @testset "Shamir network, MILPR formulation." begin
        ext = Dict(:pipe_breakpoints=>5, :pump_breakpoints=>5)
        solution = solve_wf("../test/data/epanet/shamir.inp", MILPRWaterModel, cbc, ext=ext)
        @test solution["termination_status"] == OPTIMAL
    end

    @testset "Shamir network, multinetwork MILPR formulation." begin
        data = WaterModels.parse_file("../test/data/epanet/shamir-ts.inp")
        mn_data = WaterModels.make_multinetwork(data)
        ext = Dict(:pipe_breakpoints=>5, :pump_breakpoints=>5)
        wm = instantiate_model(mn_data, MILPRWaterModel, WaterModels.build_mn_wf, ext=ext)
        solution = WaterModels.optimize_model!(wm, optimizer=cbc)
        @test solution["termination_status"] == OPTIMAL
    end

    @testset "Shamir network (with tank), multinetwork MILPR formulation." begin
        data = WaterModels.parse_file("../test/data/epanet/shamir-ts-tank.inp")
        mn_data = WaterModels.make_multinetwork(data)
        ext = Dict(:pipe_breakpoints=>5, :pump_breakpoints=>5)
        wm = instantiate_model(mn_data, MILPRWaterModel, WaterModels.build_mn_wf, ext=ext)
        solution = WaterModels.optimize_model!(wm, optimizer=cbc)
        @test solution["termination_status"] == OPTIMAL
    end

    @testset "2PRVs network, NLP formulation." begin
       data = WaterModels.parse_file("../test/data/epanet/2PRVs.inp")
       wm = instantiate_model(data, NLPWaterModel, WaterModels.build_wf)
       f = Juniper.register(head_loss_args(wm)..., autodiff=false)
       juniper = JuMP.optimizer_with_attributes(Juniper.Optimizer,
           "nl_solver"=>ipopt, "registered_functions"=>[f], "log_levels"=>[])
       solution = WaterModels.optimize_model!(wm, optimizer=juniper)

       @test solution["termination_status"] == LOCALLY_SOLVED
       @test isapprox(solution["solution"]["link"]["1"]["q"], 0.012618, rtol=1.0e-3)
       @test isapprox(solution["solution"]["link"]["4"]["q"], 0.0, atol=1.0e-4)
       @test isapprox(solution["solution"]["link"]["5"]["q"], 0.012618, rtol=1.0e-3)
       @test isapprox(solution["solution"]["node"]["1"]["h"], 45.720001, rtol=1.0e-3)
       @test isapprox(solution["solution"]["node"]["3"]["h"], 35.166771, rtol=1.0e-3)
    end

    @testset "Balerma network, NLP formulation." begin
        network = WaterModels.parse_file("../test/data/epanet/balerma.inp")
        modifications = WaterModels.parse_file("../test/data/json/balerma.json")
        InfrastructureModels.update_data!(network, modifications)
        solution = solve_wf(network, NLPWaterModel, ipopt)

        @test solution["termination_status"] == LOCALLY_SOLVED
        @test isapprox(solution["solution"]["link"]["1"]["q"], -0.002498, rtol=1.0e-3)
        @test isapprox(solution["solution"]["link"]["429"]["q"], 0.000864, rtol=1.0e-3)
        @test isapprox(solution["solution"]["link"]["540"]["q"], -0.008713, rtol=1.0e-3)
        @test isapprox(solution["solution"]["node"]["191"]["h"], 118.006958, rtol=1.0e-3)
        @test isapprox(solution["solution"]["node"]["381"]["h"], 86.70694, rtol=1.0e-3)
        @test isapprox(solution["solution"]["node"]["127001"]["h"], 84.481102, rtol=1.0e-3)
    end

    @testset "Example 1 network (with tanks and pump), NLP formulation." begin
        data = WaterModels.parse_file("../test/data/epanet/example_1-sp.inp")
        wm = instantiate_model(data, NLPWaterModel, WaterModels.build_wf)
        f = Juniper.register(head_loss_args(wm)..., autodiff=false)
        juniper = JuMP.optimizer_with_attributes(Juniper.Optimizer,
            "nl_solver"=>ipopt, "registered_functions"=>[f], "log_levels"=>[])
        solution = WaterModels.optimize_model!(wm, optimizer=juniper)
        @test solution["termination_status"] == LOCALLY_SOLVED
    end

    @testset "Example 1 network (with tanks and links), multinetwork NLP formulation." begin
        data = WaterModels.parse_file("../test/data/epanet/example_1.inp")
        mn_data = WaterModels.make_multinetwork(data)
        wm = instantiate_model(mn_data, NLPWaterModel, WaterModels.build_mn_wf)
        f = Juniper.register(head_loss_args(wm)..., autodiff=false)
        juniper = JuMP.optimizer_with_attributes(Juniper.Optimizer,
            "nl_solver"=>ipopt, "registered_functions"=>[f], "log_levels"=>[])
        solution = WaterModels.optimize_model!(wm, optimizer=juniper)
        @test solution["termination_status"] == LOCALLY_SOLVED
    end

    @testset "Richmond (single time) network, NLP formulation." begin
        data = WaterModels.parse_file("../test/data/epanet/richmond-skeleton-sp.inp")
        wm = instantiate_model(data, NLPWaterModel, WaterModels.build_wf)
        f = Juniper.register(head_loss_args(wm)..., autodiff=false)
        juniper = JuMP.optimizer_with_attributes(Juniper.Optimizer,
            "nl_solver"=>ipopt, "registered_functions"=>[f], "log_levels"=>[])
        solution = WaterModels.optimize_model!(wm, optimizer=juniper)
        @test solution["termination_status"] == LOCALLY_SOLVED
    end

    @testset "Shamir network, NLP formulation." begin
        solution = solve_wf("../test/data/epanet/shamir.inp", NLPWaterModel, ipopt)
        @test solution["termination_status"] == LOCALLY_SOLVED
        @test isapprox(solution["solution"]["link"]["2"]["q"], 0.093565, rtol=1.0e-3)
        @test isapprox(solution["solution"]["link"]["6"]["q"], 0.055710, rtol=1.0e-3)
        @test isapprox(solution["solution"]["node"]["2"]["h"], 203.247650, rtol=1.0e-3)
        @test isapprox(solution["solution"]["node"]["6"]["h"], 195.445953, rtol=1.0e-3)
        @test isapprox(solution["solution"]["reservoir"]["1"]["qr"], 0.31109, rtol=1.0e-3)
    end

    @testset "Shamir network, multinetwork NLP formulation." begin
        data = WaterModels.parse_file("../test/data/epanet/shamir-ts.inp")
        mn_data = WaterModels.make_multinetwork(data)
        wm = instantiate_model(mn_data, NLPWaterModel, WaterModels.build_mn_wf)
        solution = _IM.optimize_model!(wm, optimizer=ipopt)

        @test solution["termination_status"] == LOCALLY_SOLVED
        @test isapprox(solution["solution"]["nw"]["1"]["link"]["2"]["q"], 0.09356500, rtol=1.0e-3)
        @test isapprox(solution["solution"]["nw"]["1"]["link"]["6"]["q"], 0.05571000, rtol=1.0e-3)
        @test isapprox(solution["solution"]["nw"]["1"]["reservoir"]["1"]["qr"], 0.31109, rtol=1.0e-3)
        @test isapprox(solution["solution"]["nw"]["2"]["link"]["2"]["q"], 0.04678500, rtol=1.0e-3)
        @test isapprox(solution["solution"]["nw"]["2"]["link"]["6"]["q"], 0.02785300, rtol=1.0e-3)
        @test isapprox(solution["solution"]["nw"]["2"]["reservoir"]["1"]["qr"], 0.5*0.31109, rtol=1.0e-3)
        @test isapprox(solution["solution"]["nw"]["3"]["link"]["2"]["q"], 0.02339300, rtol=1.0e-3)
        @test isapprox(solution["solution"]["nw"]["3"]["link"]["6"]["q"], 0.01392600, rtol=1.0e-3)
        @test isapprox(solution["solution"]["nw"]["3"]["reservoir"]["1"]["qr"], 0.25*0.31109, rtol=1.0e-3)
    end

    @testset "Shamir network (with tank), multinetwork NLP formulation." begin
        data = WaterModels.parse_file("../test/data/epanet/shamir-ts-tank.inp")
        mn_data = WaterModels.make_multinetwork(data)
        wm = instantiate_model(mn_data, NLPWaterModel, WaterModels.build_mn_wf)
        f = Juniper.register(head_loss_args(wm)..., autodiff=false)
        juniper = JuMP.optimizer_with_attributes(Juniper.Optimizer,
            "nl_solver"=>ipopt, "registered_functions"=>[f], "log_levels"=>[])
        solution = _IM.optimize_model!(wm, optimizer=juniper)

        @test solution["termination_status"] == LOCALLY_SOLVED
        @test isapprox(solution["solution"]["nw"]["1"]["link"]["2"]["q"], 0.09356500, rtol=1.0e-3)
        @test isapprox(solution["solution"]["nw"]["1"]["link"]["6"]["q"], 0.05571000, rtol=1.0e-3)
        @test isapprox(solution["solution"]["nw"]["2"]["link"]["2"]["q"], 0.04678500, rtol=1.0e-3)
        @test isapprox(solution["solution"]["nw"]["2"]["link"]["6"]["q"], 0.02785300, rtol=1.0e-3)
        @test isapprox(solution["solution"]["nw"]["3"]["link"]["2"]["q"], 0.02339300, rtol=1.0e-3)
        @test isapprox(solution["solution"]["nw"]["3"]["link"]["6"]["q"], 0.01392600, rtol=1.0e-3)
        @test isapprox(solution["solution"]["nw"]["1"]["node"]["1"]["h"], 220.000000, rtol=1.0e-3)
        @test isapprox(solution["solution"]["nw"]["1"]["node"]["3"]["h"], 200.466599, rtol=1.0e-3)
        @test isapprox(solution["solution"]["nw"]["1"]["node"]["5"]["h"], 193.808380, rtol=1.0e-3)
        @test isapprox(solution["solution"]["nw"]["2"]["node"]["1"]["h"], 216.435196, rtol=1.0e-3)
        @test isapprox(solution["solution"]["nw"]["2"]["node"]["3"]["h"], 211.023941, rtol=1.0e-3)
        @test isapprox(solution["solution"]["nw"]["2"]["node"]["5"]["h"], 209.179306, rtol=1.0e-3)
        @test isapprox(solution["solution"]["nw"]["3"]["node"]["1"]["h"], 214.652771, rtol=1.0e-3)
        @test isapprox(solution["solution"]["nw"]["3"]["node"]["3"]["h"], 213.153824, rtol=1.0e-3)
        @test isapprox(solution["solution"]["nw"]["3"]["node"]["5"]["h"], 212.642853, rtol=1.0e-3)
    end
end
