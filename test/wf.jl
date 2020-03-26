@testset "Water Flow Problems" begin
    @testset "Balerma network, CNLP formulation." begin
        network = WaterModels.parse_file("../test/data/epanet/balerma.inp")
        modifications = WaterModels.parse_file("../test/data/json/balerma.json")
        InfrastructureModels.update_data!(network, modifications)
        solution = run_wf(network, CNLPWaterModel, ipopt,
            setting=Dict("output"=>Dict("duals"=>true)),
            solution_processors=[sol_data_model!])

        @test solution["termination_status"] == LOCALLY_SOLVED
        @test isapprox(solution["solution"]["link"]["1"]["q"], -0.002498, rtol=1.0e-3)
        @test isapprox(solution["solution"]["link"]["429"]["q"], 0.000864, rtol=1.0e-3)
        @test isapprox(solution["solution"]["link"]["540"]["q"], -0.008713, rtol=1.0e-3)
        @test isapprox(solution["solution"]["node"]["191"]["h"], 118.006958, rtol=1.0e-3)
        @test isapprox(solution["solution"]["node"]["381"]["h"], 86.70694, rtol=1.0e-3)
        @test isapprox(solution["solution"]["node"]["127001"]["h"], 84.481102, rtol=1.0e-3)
    end

    @testset "Example 1 network (with tanks and pumps), CNLP formulation." begin
        # The CNLP formulation does not support tanks, so an error is thrown.
        @test_throws ErrorException run_wf("../test/data/epanet/example_1-sp.inp", CNLPWaterModel, ipopt)
    end

    @testset "Example 1 network (with tanks and pumps), multinetwork CNLP formulation." begin
        data = WaterModels.parse_file("../test/data/epanet/example_1.inp")
        mn_data = WaterModels.make_multinetwork(data)

        # The CNLP formulation does not support tanks, so an error is thrown.
        @test_throws ErrorException instantiate_model(mn_data, CNLPWaterModel, WaterModels.post_mn_wf)
    end

    @testset "Richmond (single time) network, CNLP formulation." begin
        # The CNLP formulation does not support tanks, so an error is thrown.
        @test_throws ErrorException run_wf("../test/data/epanet/richmond-skeleton-sp.inp", CNLPWaterModel, ipopt)
    end

    @testset "Richmond network, multinetwork CNLP formulation." begin
        data = WaterModels.parse_file("../test/data/epanet/richmond-skeleton.inp")
        mn_data = WaterModels.make_multinetwork(data)

        # The CNLP formulation does not support tanks, so an error is thrown.
        @test_throws ErrorException instantiate_model(mn_data, CNLPWaterModel, WaterModels.post_mn_wf)
    end

    @testset "Shamir network, CNLP formulation." begin
        solution = run_wf("../test/data/epanet/shamir.inp", CNLPWaterModel,
            ipopt, setting=Dict("output"=>Dict("duals"=>true)),
            solution_processors=[sol_data_model!])
        @test solution["termination_status"] == LOCALLY_SOLVED
        @test isapprox(solution["solution"]["link"]["2"]["q"], 0.093565, rtol=1.0e-3)
        @test isapprox(solution["solution"]["link"]["6"]["q"], 0.055710, rtol=1.0e-3)
        @test isapprox(solution["solution"]["node"]["2"]["h"], 203.247650, rtol=1.0e-3)
        @test isapprox(solution["solution"]["node"]["6"]["h"], 195.445953, rtol=1.0e-3)
        @test isapprox(solution["solution"]["reservoir"]["1"]["qr"], 0.31109, rtol=1.0e-3)
    end

    @testset "Shamir network, multinetwork CNLP formulation." begin
        data = WaterModels.parse_file("../test/data/epanet/shamir-ts.inp")
        mn_data = WaterModels.make_multinetwork(data)
        wm = instantiate_model(mn_data, CNLPWaterModel, WaterModels.post_mn_wf,
            setting=Dict("output"=>Dict("duals"=>true)))
        solution = _IM.optimize_model!(wm, optimizer=ipopt,
            solution_processors=[sol_data_model!])

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

    @testset "Shamir network (with tank), multinetwork CNLP formulation." begin
        data = WaterModels.parse_file("../test/data/epanet/shamir-ts-tank.inp")
        mn_data = WaterModels.make_multinetwork(data)

        # The CNLP formulation does not support tanks, so an error is thrown.
        @test_throws ErrorException instantiate_model(mn_data, CNLPWaterModel, WaterModels.post_mn_wf)
    end

    @testset "Balerma network, MICP formulation." begin
        data = WaterModels.parse_file("../test/data/epanet/balerma.inp")
        modifications = WaterModels.parse_file("../test/data/json/balerma.json")
        InfrastructureModels.update_data!(data, modifications)
        wm = instantiate_model(data, MICPWaterModel, WaterModels.post_wf)
        f = Juniper.register(head_loss_args(wm)..., autodiff=false)
        juniper = JuMP.optimizer_with_attributes(Juniper.Optimizer, "nl_solver"=>ipopt,
            "registered_functions"=>[f], "log_levels"=>[])

        # TODO: The below takes too long to execute.
        #solution = WaterModels.optimize_model!(wm, optimizer=juniper)
        #@test solution["termination_status"] == LOCALLY_SOLVED
    end

    @testset "Example 1 network (with tanks and pumps), MICP formulation." begin
        data = WaterModels.parse_file("../test/data/epanet/example_1-sp.inp")
        wm = instantiate_model(data, MICPWaterModel, WaterModels.post_wf)
        f = Juniper.register(head_loss_args(wm)..., autodiff=false)
        juniper = JuMP.optimizer_with_attributes(Juniper.Optimizer,
            "nl_solver"=>ipopt, "registered_functions"=>[f], "log_levels"=>[])
        solution = WaterModels.optimize_model!(wm, optimizer=juniper)
        @test solution["termination_status"] == LOCALLY_SOLVED
    end

    @testset "Example 1 network (with tanks and pumps), multinetwork MICP formulation." begin
        data = WaterModels.parse_file("../test/data/epanet/example_1.inp")
        mn_data = WaterModels.make_multinetwork(data)
        wm = instantiate_model(mn_data, MICPWaterModel, WaterModels.post_mn_wf)
        f = Juniper.register(head_loss_args(wm)..., autodiff=false)
        juniper = JuMP.optimizer_with_attributes(Juniper.Optimizer,
            "nl_solver"=>ipopt, "registered_functions"=>[f], "log_levels"=>[])
        solution = WaterModels.optimize_model!(wm, optimizer=juniper)
        @test solution["termination_status"] == LOCALLY_SOLVED
    end

    @testset "Richmond (single time) network, MICP formulation." begin
        data = WaterModels.parse_file("../test/data/epanet/richmond-skeleton-sp.inp")
        wm = instantiate_model(data, MICPWaterModel, WaterModels.post_wf)
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
        wm = instantiate_model(mn_data, MICPWaterModel, WaterModels.post_mn_wf)
        f = Juniper.register(head_loss_args(wm)..., autodiff=false)
        juniper = JuMP.optimizer_with_attributes(Juniper.Optimizer,
            "nl_solver"=>ipopt, "registered_functions"=>[f], "log_levels"=>[])
        # TODO: The below takes too long to execute.
        #solution = WaterModels.optimize_model!(wm, optimizer=juniper)
        #@test solution["termination_status"] == LOCALLY_SOLVED
    end

    @testset "Shamir network, MICP formulation." begin
        data = WaterModels.parse_file("../test/data/epanet/shamir.inp")
        wm = instantiate_model(data, MICPWaterModel, WaterModels.post_wf)
        f = Juniper.register(head_loss_args(wm)..., autodiff=false)
        juniper = JuMP.optimizer_with_attributes(Juniper.Optimizer,
            "nl_solver"=>ipopt, "registered_functions"=>[f], "log_levels"=>[])
        solution = WaterModels.optimize_model!(wm, optimizer=juniper)
        @test solution["termination_status"] == LOCALLY_SOLVED
    end

    @testset "Shamir network, multinetwork MICP formulation." begin
        data = WaterModels.parse_file("../test/data/epanet/shamir-ts.inp")
        mn_data = WaterModels.make_multinetwork(data)
        wm = instantiate_model(mn_data, MICPWaterModel, WaterModels.post_mn_wf)
        f = Juniper.register(head_loss_args(wm)..., autodiff=false)
        juniper = JuMP.optimizer_with_attributes(Juniper.Optimizer,
            "nl_solver"=>ipopt, "registered_functions"=>[f], "log_levels"=>[])
        solution = WaterModels.optimize_model!(wm, optimizer=juniper)
        @test solution["termination_status"] == LOCALLY_SOLVED
    end

    @testset "Shamir network (with tank), multinetwork MICP formulation." begin
        data = WaterModels.parse_file("../test/data/epanet/shamir-ts-tank.inp")
        mn_data = WaterModels.make_multinetwork(data)
        wm = instantiate_model(mn_data, MICPWaterModel, WaterModels.post_mn_wf)
        f = Juniper.register(head_loss_args(wm)..., autodiff=false)
        juniper = JuMP.optimizer_with_attributes(Juniper.Optimizer,
            "nl_solver"=>ipopt, "registered_functions"=>[f], "log_levels"=>[])
        solution = WaterModels.optimize_model!(wm, optimizer=juniper)
        @test solution["termination_status"] == LOCALLY_SOLVED
    end

    @testset "Balerma network, MILPR formulation." begin
        data = WaterModels.parse_file("../test/data/epanet/balerma.inp")
        modifications = WaterModels.parse_file("../test/data/json/balerma.json")
        InfrastructureModels.update_data!(data, modifications)
        ext = Dict(:num_breakpoints => 5)
        wm = instantiate_model(data, MILPRWaterModel, WaterModels.post_wf, ext=ext)
        solution = WaterModels.optimize_model!(wm, optimizer=cbc)
        @test solution["termination_status"] == OPTIMAL
    end

    @testset "Example 1 network (with tanks and pumps), MILPR formulation." begin
        ext = Dict(:num_breakpoints => 5)
        data = WaterModels.parse_file("../test/data/epanet/example_1-sp.inp")
        wm = instantiate_model(data, MILPRWaterModel, WaterModels.post_wf, ext=ext)
        juniper = JuMP.optimizer_with_attributes(Juniper.Optimizer, "nl_solver"=>ipopt, "log_levels"=>[])
        solution = WaterModels.optimize_model!(wm, optimizer=juniper)
        @test solution["termination_status"] == LOCALLY_SOLVED
    end

    @testset "Example 1 network (with tanks and pumps), multinetwork MILPR formulation." begin
        ext = Dict(:num_breakpoints => 5)
        data = WaterModels.parse_file("../test/data/epanet/example_1.inp")
        mn_data = WaterModels.make_multinetwork(data)
        wm = instantiate_model(mn_data, MILPRWaterModel, WaterModels.post_mn_wf, ext=ext)
        juniper = JuMP.optimizer_with_attributes(Juniper.Optimizer, "nl_solver"=>ipopt, "log_levels"=>[])
        solution = WaterModels.optimize_model!(wm, optimizer=juniper)
        @test solution["termination_status"] == LOCALLY_SOLVED
    end

    @testset "Richmond (single time) network, MILPR formulation." begin
        ext = Dict(:num_breakpoints => 5)
        data = WaterModels.parse_file("../test/data/epanet/richmond-skeleton-sp.inp")
        wm = instantiate_model(data, MILPRWaterModel, WaterModels.post_wf, ext=ext)
        juniper = JuMP.optimizer_with_attributes(Juniper.Optimizer, "nl_solver"=>ipopt, "log_levels"=>[])
        # TODO: The below takes too long to execute.
        #solution = WaterModels.optimize_model!(wm, optimizer=juniper)
        #@test solution["termination_status"] == LOCALLY_SOLVED
    end

    @testset "Richmond network, multinetwork MILPR formulation." begin
        ext = Dict(:num_breakpoints => 5)
        data = WaterModels.parse_file("../test/data/epanet/richmond-skeleton.inp")
        mn_data = WaterModels.make_multinetwork(data)
        wm = instantiate_model(mn_data, MILPRWaterModel, WaterModels.post_mn_wf, ext=ext)
        juniper = JuMP.optimizer_with_attributes(Juniper.Optimizer, "nl_solver"=>ipopt, "log_levels"=>[])
        # TODO: The below takes too long to execute.
        #solution = WaterModels.optimize_model!(wm, optimizer=juniper)
        #@test solution["termination_status"] == LOCALLY_SOLVED
    end

    @testset "Shamir network, MILPR formulation." begin
        ext = Dict(:num_breakpoints => 5)
        solution = run_wf("../test/data/epanet/shamir.inp", MILPRWaterModel, cbc, ext=ext)
        @test solution["termination_status"] == OPTIMAL
    end

    @testset "Shamir network, multinetwork MILPR formulation." begin
        data = WaterModels.parse_file("../test/data/epanet/shamir-ts.inp")
        mn_data = WaterModels.make_multinetwork(data)
        ext = Dict(:num_breakpoints => 5)
        wm = instantiate_model(mn_data, MILPRWaterModel, WaterModels.post_mn_wf, ext=ext)
        solution = WaterModels.optimize_model!(wm, optimizer=cbc)
        @test solution["termination_status"] == OPTIMAL
    end

    @testset "Shamir network (with tank), multinetwork MILPR formulation." begin
        data = WaterModels.parse_file("../test/data/epanet/shamir-ts-tank.inp")
        mn_data = WaterModels.make_multinetwork(data)
        ext = Dict(:num_breakpoints => 5)
        wm = instantiate_model(mn_data, MILPRWaterModel, WaterModels.post_mn_wf, ext=ext)
        solution = WaterModels.optimize_model!(wm, optimizer=cbc)
        @test solution["termination_status"] == OPTIMAL
    end

    @testset "Balerma network, NCNLP formulation." begin
        network = WaterModels.parse_file("../test/data/epanet/balerma.inp")
        modifications = WaterModels.parse_file("../test/data/json/balerma.json")
        InfrastructureModels.update_data!(network, modifications)
        solution = run_wf(network, NCNLPWaterModel, ipopt)

        @test solution["termination_status"] == LOCALLY_SOLVED
        @test isapprox(solution["solution"]["link"]["1"]["q"], -0.002498, rtol=1.0e-3)
        @test isapprox(solution["solution"]["link"]["429"]["q"], 0.000864, rtol=1.0e-3)
        @test isapprox(solution["solution"]["link"]["540"]["q"], -0.008713, rtol=1.0e-3)
        @test isapprox(solution["solution"]["node"]["191"]["h"], 118.006958, rtol=1.0e-3)
        @test isapprox(solution["solution"]["node"]["381"]["h"], 86.70694, rtol=1.0e-3)
        @test isapprox(solution["solution"]["node"]["127001"]["h"], 84.481102, rtol=1.0e-3)
    end

    @testset "Example 1 network (with tanks and pump), NCNLP formulation." begin
        data = WaterModels.parse_file("../test/data/epanet/example_1-sp.inp")
        wm = instantiate_model(data, NCNLPWaterModel, WaterModels.post_wf)
        f = Juniper.register(head_loss_args(wm)..., autodiff=false)
        juniper = JuMP.optimizer_with_attributes(Juniper.Optimizer,
            "nl_solver"=>ipopt, "registered_functions"=>[f], "log_levels"=>[])
        solution = WaterModels.optimize_model!(wm, optimizer=juniper)

        @test solution["termination_status"] == LOCALLY_SOLVED
        @test isapprox(solution["solution"]["link"]["9"]["q"], 0.117737, rtol=1.0e-3)
        @test isapprox(solution["solution"]["link"]["110"]["q"], -0.048338, rtol=1.0e-3)
        @test isapprox(solution["solution"]["link"]["122"]["q"], 0.003734, rtol=1.0e-3)
        @test isapprox(solution["solution"]["node"]["9"]["h"], 243.839996, rtol=1.0e-3)
        @test isapprox(solution["solution"]["node"]["10"]["h"], 306.125092, rtol=1.0e-3)
        @test isapprox(solution["solution"]["node"]["23"]["h"], 295.243073, rtol=1.0e-3)
    end

    @testset "Example 1 network (with tanks and links), multinetwork NCNLP formulation." begin
        data = WaterModels.parse_file("../test/data/epanet/example_1.inp")
        mn_data = WaterModels.make_multinetwork(data)
        wm = instantiate_model(mn_data, NCNLPWaterModel, WaterModels.post_mn_wf)
        f = Juniper.register(head_loss_args(wm)..., autodiff=false)
        juniper = JuMP.optimizer_with_attributes(Juniper.Optimizer,
            "nl_solver"=>ipopt, "registered_functions"=>[f], "log_levels"=>[])
        solution = WaterModels.optimize_model!(wm, optimizer=juniper)

        @test solution["termination_status"] == LOCALLY_SOLVED
        @test isapprox(solution["solution"]["nw"]["1"]["link"]["9"]["q"], 0.117737, rtol=1.0e-3)
        @test isapprox(solution["solution"]["nw"]["1"]["link"]["110"]["q"], -0.048338, rtol=1.0e-3)
        @test isapprox(solution["solution"]["nw"]["1"]["link"]["122"]["q"], 0.003734, rtol=1.0e-3)
        @test isapprox(solution["solution"]["nw"]["1"]["node"]["9"]["h"], 243.839996, rtol=1.0e-3)
        @test isapprox(solution["solution"]["nw"]["1"]["node"]["10"]["h"], 306.125092, rtol=1.0e-3)
        @test isapprox(solution["solution"]["nw"]["1"]["node"]["23"]["h"], 295.243073, rtol=1.0e-3)
        @test isapprox(solution["solution"]["nw"]["3"]["link"]["9"]["q"], 0.115926, rtol=1.0e-3)
        @test isapprox(solution["solution"]["nw"]["3"]["link"]["110"]["q"], -0.032647, rtol=1.0e-3)
        @test isapprox(solution["solution"]["nw"]["3"]["link"]["122"]["q"], 0.004907, rtol=1.0e-3)
        @test isapprox(solution["solution"]["nw"]["3"]["node"]["9"]["h"], 243.839996, rtol=1.0e-3)
        @test isapprox(solution["solution"]["nw"]["3"]["node"]["10"]["h"], 307.325653, rtol=1.0e-3)
        @test isapprox(solution["solution"]["nw"]["3"]["node"]["23"]["h"], 296.683014, rtol=1.0e-3)
    end

    @testset "Richmond (single time) network, NCNLP formulation." begin
        data = WaterModels.parse_file("../test/data/epanet/richmond-skeleton-sp.inp")
        wm = instantiate_model(data, NCNLPWaterModel, WaterModels.post_wf)
        f = Juniper.register(head_loss_args(wm)..., autodiff=false)
        juniper = JuMP.optimizer_with_attributes(Juniper.Optimizer,
            "nl_solver"=>ipopt, "registered_functions"=>[f], "log_levels"=>[])
        solution = WaterModels.optimize_model!(wm, optimizer=juniper)

        @test solution["termination_status"] == LOCALLY_SOLVED
        @test isapprox(solution["solution"]["node"]["17"]["h"], 243.052536, rtol=1.0e-3)
        @test isapprox(solution["solution"]["node"]["39"]["h"], 187.250000, rtol=1.0e-3)
        @test isapprox(solution["solution"]["node"]["45"]["h"], 243.120010, rtol=1.0e-3)
        @test isapprox(solution["solution"]["link"]["18"]["q"], 5.308596e-03, rtol=1.0e-3)
        @test isapprox(solution["solution"]["link"]["44"]["q"], 9.160000e-03, rtol=1.0e-3)
        @test isapprox(solution["solution"]["link"]["47"]["q"], 0.00000000, atol=1.0e-7)
    end

    @testset "Richmond network, multinetwork NCNLP formulation." begin
        data = WaterModels.parse_file("../test/data/epanet/richmond-skeleton.inp")
        mn_data = WaterModels.make_multinetwork(data)
        wm = instantiate_model(mn_data, NCNLPWaterModel, WaterModels.post_mn_wf)
        f = Juniper.register(head_loss_args(wm)..., autodiff=false)
        juniper = JuMP.optimizer_with_attributes(Juniper.Optimizer,
            "nl_solver"=>ipopt, "registered_functions"=>[f], "log_levels"=>[])
        solution = WaterModels.optimize_model!(wm, optimizer=juniper)

        # TODO: The below does not solve using Juniper.
        # @test solution["termination_status"] == LOCALLY_SOLVED
        # @test isapprox(solution["solution"]["nw"]["1"]["node"]["17"]["h"], 243.052536, rtol=1.0e-3)
        # @test isapprox(solution["solution"]["nw"]["1"]["node"]["39"]["h"], 187.250000, rtol=1.0e-3)
        # @test isapprox(solution["solution"]["nw"]["1"]["node"]["45"]["h"], 243.120010, rtol=1.0e-3)
        # @test isapprox(solution["solution"]["nw"]["1"]["link"]["18"]["q"], 5.308596e-03, rtol=1.0e-3)
        # @test isapprox(solution["solution"]["nw"]["1"]["link"]["44"]["q"], 9.160000e-03, rtol=1.0e-3)
        # @test isapprox(solution["solution"]["nw"]["1"]["link"]["47"]["q"], 0.00000000, atol=1.0e-7)
    end

    @testset "Shamir network, NCNLP formulation." begin
        solution = run_wf("../test/data/epanet/shamir.inp", NCNLPWaterModel, ipopt)
        @test solution["termination_status"] == LOCALLY_SOLVED
        @test isapprox(solution["solution"]["link"]["2"]["q"], 0.093565, rtol=1.0e-3)
        @test isapprox(solution["solution"]["link"]["6"]["q"], 0.055710, rtol=1.0e-3)
        @test isapprox(solution["solution"]["node"]["2"]["h"], 203.247650, rtol=1.0e-3)
        @test isapprox(solution["solution"]["node"]["6"]["h"], 195.445953, rtol=1.0e-3)
        @test isapprox(solution["solution"]["reservoir"]["1"]["qr"], 0.31109, rtol=1.0e-3)
    end

    @testset "Shamir network, multinetwork NCNLP formulation." begin
        data = WaterModels.parse_file("../test/data/epanet/shamir-ts.inp")
        mn_data = WaterModels.make_multinetwork(data)
        wm = instantiate_model(mn_data, NCNLPWaterModel, WaterModels.post_mn_wf)
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

    @testset "Shamir network (with tank), multinetwork NCNLP formulation." begin
        data = WaterModels.parse_file("../test/data/epanet/shamir-ts-tank.inp")
        mn_data = WaterModels.make_multinetwork(data)
        wm = instantiate_model(mn_data, NCNLPWaterModel, WaterModels.post_mn_wf)
        solution = _IM.optimize_model!(wm, optimizer=ipopt)

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
