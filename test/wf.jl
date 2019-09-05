@testset "Water Flow Problems" begin
    # Set test network paths.
    balerma_path = "../test/data/epanet/balerma.inp"
    example_1_path = "../test/data/epanet/example_1.inp"
    example_1_sp_path = "../test/data/epanet/example_1-sp.inp"
    richmond_skeleton_path = "../test/data/epanet/richmond-skeleton.inp"
    richmond_skeleton_sp_path = "../test/data/epanet/richmond-skeleton-sp.inp"
    shamir_path = "../test/data/epanet/shamir.inp"
    shamir_ts_path = "../test/data/epanet/shamir-ts.inp"
    shamir_ts_tank_path = "../test/data/epanet/shamir-ts-tank.inp"

    @testset "Balerma network (unknown flow directions), CNLP formulation." begin
        solution = run_wf(balerma_path, CNLPWaterModel, ipopt)
        @test solution["termination_status"] == LOCALLY_SOLVED
    end

    #@testset "Richmond network, multinetwork NCNLP formulation." begin
    #    network_data = WaterModels.parse_file(richmond_skeleton_path)
    #    mn_data = WaterModels.make_multinetwork(network_data)
    #    wm = build_model(mn_data, NCNLPWaterModel, WaterModels.post_mn_wf, multinetwork=true)
    #    f_1, f_2, f_3, f_4, f_5 = fun(wm, :head_loss)
    #    f = Juniper.register(f_1, f_2, f_3, f_4, f_5, autodiff=false)
    #    juniper = JuMP.with_optimizer(Juniper.Optimizer, nl_solver=ipopt, registered_functions=[f]) #, log_levels=[])
    #    solution = WaterModels.optimize_model!(wm, juniper)

    #    @test solution["termination_status"] == LOCALLY_SOLVED
    #    @test isapprox(solution["solution"]["nw"]["1"]["node"]["17"]["h"], 243.052536, rtol=1.0e-3)
    #    @test isapprox(solution["solution"]["nw"]["1"]["node"]["39"]["h"], 187.250000, rtol=1.0e-3)
    #    @test isapprox(solution["solution"]["nw"]["1"]["node"]["45"]["h"], 243.120010, rtol=1.0e-3)
    #    @test isapprox(solution["solution"]["nw"]["1"]["pipe"]["18"]["q"], 5.308596e-03, rtol=1.0e-3)
    #    @test isapprox(solution["solution"]["nw"]["1"]["pipe"]["44"]["q"], 9.160000e-03, rtol=1.0e-3)
    #    @test isapprox(solution["solution"]["nw"]["1"]["pump"]["47"]["q"], 0.00000000, atol=1.0e-7)
    #end

    @testset "Single-time Richmond network, NCNLP formulation." begin
        network_data = WaterModels.parse_file(richmond_skeleton_sp_path)
        wm = build_model(network_data, NCNLPWaterModel, WaterModels.post_wf, multinetwork=false)
        f = Juniper.register(fun(wm, :head_loss)..., autodiff=false)
        juniper = JuMP.with_optimizer(Juniper.Optimizer, nl_solver=ipopt, registered_functions=[f], log_levels=[])
        solution = WaterModels.optimize_model!(wm, juniper)

        @test solution["termination_status"] == LOCALLY_SOLVED
        @test isapprox(solution["solution"]["node"]["17"]["h"], 243.052536, rtol=1.0e-3)
        @test isapprox(solution["solution"]["node"]["39"]["h"], 187.250000, rtol=1.0e-3)
        @test isapprox(solution["solution"]["node"]["45"]["h"], 243.120010, rtol=1.0e-3)
        @test isapprox(solution["solution"]["pipe"]["18"]["q"], 5.308596e-03, rtol=1.0e-3)
        @test isapprox(solution["solution"]["pipe"]["44"]["q"], 9.160000e-03, rtol=1.0e-3)
        @test isapprox(solution["solution"]["pump"]["47"]["q"], 0.00000000, atol=1.0e-7)
    end

    @testset "Shamir network (unknown flow directions), CNLP formulation." begin
        solution = run_wf(shamir_path, CNLPWaterModel, ipopt)
        @test solution["termination_status"] == LOCALLY_SOLVED
        @test isapprox(solution["solution"]["pipe"]["2"]["q"], 0.093565, rtol=1.0e-3)
        @test isapprox(solution["solution"]["pipe"]["6"]["q"], 0.055710, rtol=1.0e-3)
    end

    @testset "Shamir network, CNLP formulation." begin
        solution = run_wf(shamir_path, CNLPWaterModel, ipopt)
        @test solution["termination_status"] == LOCALLY_SOLVED
        @test isapprox(solution["solution"]["pipe"]["2"]["q"], 0.093565, rtol=1.0e-3)
        @test isapprox(solution["solution"]["pipe"]["6"]["q"], 0.055710, rtol=1.0e-3)
        @test isapprox(solution["solution"]["node"]["2"]["h"], 203.247650, rtol=1.0e-3)
        @test isapprox(solution["solution"]["node"]["6"]["h"], 195.445953, rtol=1.0e-3)
    end

    @testset "Shamir network (unknown flow directions), NCNLP formulation." begin
        solution = run_wf(shamir_path, NCNLPWaterModel, ipopt)
        @test solution["termination_status"] == LOCALLY_SOLVED
        @test isapprox(solution["solution"]["pipe"]["2"]["q"], 0.093565, rtol=1.0e-3)
        @test isapprox(solution["solution"]["pipe"]["6"]["q"], 0.055710, rtol=1.0e-3)
        @test isapprox(solution["solution"]["node"]["2"]["h"], 203.247650, rtol=1.0e-3)
        @test isapprox(solution["solution"]["node"]["6"]["h"], 195.445953, rtol=1.0e-3)
    end

    @testset "Shamir network, multinetwork NCNLP formulation." begin
        network_data = WaterModels.parse_file(shamir_ts_path)
        mn_data = WaterModels.make_multinetwork(network_data)
        wm = build_model(mn_data, NCNLPWaterModel, WaterModels.post_mn_wf, multinetwork=true)
        solution = optimize_model!(wm, ipopt)

        @test solution["termination_status"] == LOCALLY_SOLVED
        @test isapprox(solution["solution"]["nw"]["1"]["pipe"]["2"]["q"], 0.09356500, rtol=1.0e-3)
        @test isapprox(solution["solution"]["nw"]["1"]["pipe"]["6"]["q"], 0.05571000, rtol=1.0e-3)
        @test isapprox(solution["solution"]["nw"]["2"]["pipe"]["2"]["q"], 0.04678500, rtol=1.0e-3)
        @test isapprox(solution["solution"]["nw"]["2"]["pipe"]["6"]["q"], 0.02785300, rtol=1.0e-3)
        @test isapprox(solution["solution"]["nw"]["3"]["pipe"]["2"]["q"], 0.02339300, rtol=1.0e-3)
        @test isapprox(solution["solution"]["nw"]["3"]["pipe"]["6"]["q"], 0.01392600, rtol=1.0e-3)
    end

    @testset "Shamir network (with tank), multinetwork NCNLP formulation." begin
        network_data = WaterModels.parse_file(shamir_ts_tank_path)
        mn_data = WaterModels.make_multinetwork(network_data)
        wm = build_model(mn_data, NCNLPWaterModel, WaterModels.post_mn_wf, multinetwork=true)
        solution = optimize_model!(wm, ipopt)

        @test solution["termination_status"] == LOCALLY_SOLVED
        @test isapprox(solution["solution"]["nw"]["1"]["pipe"]["2"]["q"], 0.09356500, rtol=1.0e-3)
        @test isapprox(solution["solution"]["nw"]["1"]["pipe"]["6"]["q"], 0.05571000, rtol=1.0e-3)
        @test isapprox(solution["solution"]["nw"]["2"]["pipe"]["2"]["q"], 0.04678500, rtol=1.0e-3)
        @test isapprox(solution["solution"]["nw"]["2"]["pipe"]["6"]["q"], 0.02785300, rtol=1.0e-3)
        @test isapprox(solution["solution"]["nw"]["3"]["pipe"]["2"]["q"], 0.02339300, rtol=1.0e-3)
        @test isapprox(solution["solution"]["nw"]["3"]["pipe"]["6"]["q"], 0.01392600, rtol=1.0e-3)
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

    @testset "Example 1 network (with tanks and pumps), NCNLP formulation." begin
        network_data = WaterModels.parse_file(example_1_sp_path)
        wm = build_model(network_data, NCNLPWaterModel, WaterModels.post_wf)
        f = Juniper.register(fun(wm, :head_loss)..., autodiff=false)
        juniper = JuMP.with_optimizer(Juniper.Optimizer, nl_solver=ipopt, registered_functions=[f], log_levels=[])
        solution = WaterModels.optimize_model!(wm, juniper)

        @test solution["termination_status"] == LOCALLY_SOLVED
        @test isapprox(solution["solution"]["pump"]["9"]["q"], 0.117737, rtol=1.0e-3)
        @test isapprox(solution["solution"]["pipe"]["110"]["q"], -0.048338, rtol=1.0e-3)
        @test isapprox(solution["solution"]["pipe"]["122"]["q"], 0.003734, rtol=1.0e-3)
        @test isapprox(solution["solution"]["node"]["9"]["h"], 243.839996, rtol=1.0e-3)
        @test isapprox(solution["solution"]["node"]["10"]["h"], 306.125092, rtol=1.0e-3)
        @test isapprox(solution["solution"]["node"]["23"]["h"], 295.243073, rtol=1.0e-3)
    end

    @testset "Example 1 network (with tanks and pumps), multinetwork NCNLP formulation." begin
        network_data = WaterModels.parse_file(example_1_path)
        mn_data = WaterModels.make_multinetwork(network_data)
        wm = build_model(mn_data, NCNLPWaterModel, WaterModels.post_mn_wf, multinetwork=true)
        f = Juniper.register(fun(wm, :head_loss)..., autodiff=false)
        juniper = JuMP.with_optimizer(Juniper.Optimizer, nl_solver=ipopt, registered_functions=[f], log_levels=[])
        solution = WaterModels.optimize_model!(wm, juniper)

        @test solution["termination_status"] == LOCALLY_SOLVED
        @test isapprox(solution["solution"]["nw"]["1"]["pump"]["9"]["q"], 0.117737, rtol=1.0e-3)
        @test isapprox(solution["solution"]["nw"]["1"]["pipe"]["110"]["q"], -0.048338, rtol=1.0e-3)
        @test isapprox(solution["solution"]["nw"]["1"]["pipe"]["122"]["q"], 0.003734, rtol=1.0e-3)
        @test isapprox(solution["solution"]["nw"]["1"]["node"]["9"]["h"], 243.839996, rtol=1.0e-3)
        @test isapprox(solution["solution"]["nw"]["1"]["node"]["10"]["h"], 306.125092, rtol=1.0e-3)
        @test isapprox(solution["solution"]["nw"]["1"]["node"]["23"]["h"], 295.243073, rtol=1.0e-3)
        @test isapprox(solution["solution"]["nw"]["3"]["pump"]["9"]["q"], 0.115926, rtol=1.0e-3)
        @test isapprox(solution["solution"]["nw"]["3"]["pipe"]["110"]["q"], -0.032647, rtol=1.0e-3)
        @test isapprox(solution["solution"]["nw"]["3"]["pipe"]["122"]["q"], 0.004907, rtol=1.0e-3)
        @test isapprox(solution["solution"]["nw"]["3"]["node"]["9"]["h"], 243.839996, rtol=1.0e-3)
        @test isapprox(solution["solution"]["nw"]["3"]["node"]["10"]["h"], 307.325653, rtol=1.0e-3)
        @test isapprox(solution["solution"]["nw"]["3"]["node"]["23"]["h"], 296.683014, rtol=1.0e-3)
    end

    @testset "Shamir network (unknown flow directions), MICP formulation." begin
        network_data = WaterModels.parse_file(shamir_path)
        wm = build_model(network_data, MICPWaterModel, WaterModels.post_wf, multinetwork=false)
        f = Juniper.register(fun(wm, :head_loss)..., autodiff=false)
        juniper = JuMP.with_optimizer(Juniper.Optimizer, nl_solver=ipopt, registered_functions=[f], log_levels=[])
        solution = WaterModels.optimize_model!(wm, juniper)
        @test solution["termination_status"] == LOCALLY_SOLVED
    end

    @testset "Shamir network (unknown flow directions), MILPR formulation." begin
        solution = run_wf(shamir_path, MILPRWaterModel, cbc)
        @test solution["termination_status"] == OPTIMAL
        @test isapprox(solution["solution"]["pipe"]["2"]["q"], 0.093565, rtol=0.05)
        @test isapprox(solution["solution"]["pipe"]["6"]["q"], 0.055710, rtol=0.05)
        @test isapprox(solution["solution"]["node"]["2"]["h"], 203.247650, rtol=0.05)
        @test isapprox(solution["solution"]["node"]["6"]["h"], 195.445953, rtol=0.05)
    end
end
