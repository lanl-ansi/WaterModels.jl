@testset "Water Flow Problems" begin
    # Set test network paths.
    balerma_path = "../test/data/epanet/balerma.inp"
    richmond_skeleton_path = "../test/data/epanet/richmond-skeleton.inp"
    richmond_skeleton_sp_path = "../test/data/epanet/richmond-skeleton-sp.inp"
    shamir_path = "../test/data/epanet/shamir.inp"
    shamir_ts_path = "../test/data/epanet/shamir-ts.inp"
    shamir_ts_tank_path = "../test/data/epanet/shamir-ts-tank.inp"

    @testset "Balerma network (unknown flow directions), CNLP formulation." begin
        solution = run_wf(balerma_path, CNLPWaterModel, ipopt)
        @test solution["termination_status"] == MOI.LOCALLY_SOLVED
    end

    #@testset "Richmond network (unknown flow directions), NCNLP formulation." begin
    #    solution = run_wf(richmond_skeleton_path, NCNLPWaterModel, ipopt)
    #end

    #@testset "Single-time Richmond network, NCNLP formulation." begin
    #    solution = run_wf(richmond_skeleton_sp_path, NCNLPWaterModel, ipopt)
    #end

    @testset "Shamir network (unknown flow directions), CNLP formulation." begin
        solution = run_wf(shamir_path, CNLPWaterModel, ipopt)
        @test solution["termination_status"] == MOI.LOCALLY_SOLVED
        @test isapprox(solution["solution"]["pipes"]["2"]["q"], 0.093565, rtol=1.0e-4)
        @test isapprox(solution["solution"]["pipes"]["6"]["q"], 0.055710, rtol=1.0e-4)
    end

    @testset "Shamir network, multinetwork CNLP formulation." begin
        solution = run_wf(shamir_path, CNLPWaterModel, ipopt)
        @test solution["termination_status"] == MOI.LOCALLY_SOLVED
        @test isapprox(solution["solution"]["pipes"]["2"]["q"], 0.093565, rtol=1.0e-4)
        @test isapprox(solution["solution"]["pipes"]["6"]["q"], 0.055710, rtol=1.0e-4)
    end

    @testset "Shamir network (unknown flow directions), NCNLP formulation." begin
        solution = run_wf(shamir_path, NCNLPWaterModel, ipopt)
        @test solution["termination_status"] == MOI.LOCALLY_SOLVED
        @test isapprox(solution["solution"]["pipes"]["2"]["q"], 0.093565, rtol=1.0e-4)
        @test isapprox(solution["solution"]["pipes"]["6"]["q"], 0.055710, rtol=1.0e-4)
        @test isapprox(solution["solution"]["nodes"]["2"]["h"], 203.247650, rtol=1.0e-4)
        @test isapprox(solution["solution"]["nodes"]["6"]["h"], 195.445953, rtol=1.0e-4)
    end

    @testset "Shamir network, multinetwork NCNLP formulation." begin
        network_data = WaterModels.parse_file(shamir_ts_path)
        mn_data = WaterModels.make_multinetwork(network_data)
        wm = build_generic_model(mn_data, NCNLPWaterModel, WaterModels.post_mn_wf, multinetwork=true)
        solution = solve_generic_model(wm, ipopt)

        @test solution["termination_status"] == MOI.LOCALLY_SOLVED
        @test isapprox(solution["solution"]["nw"]["1"]["pipes"]["2"]["q"], 0.09356500, rtol=1.0e-4)
        @test isapprox(solution["solution"]["nw"]["1"]["pipes"]["6"]["q"], 0.05571000, rtol=1.0e-4)
        @test isapprox(solution["solution"]["nw"]["2"]["pipes"]["2"]["q"], 0.04678500, rtol=1.0e-4)
        @test isapprox(solution["solution"]["nw"]["2"]["pipes"]["6"]["q"], 0.02785300, rtol=1.0e-4)
        @test isapprox(solution["solution"]["nw"]["3"]["pipes"]["2"]["q"], 0.02339300, rtol=1.0e-4)
        @test isapprox(solution["solution"]["nw"]["3"]["pipes"]["6"]["q"], 0.01392600, rtol=1.0e-4)
    end

    @testset "Shamir network (with tank), multinetwork NCNLP formulation." begin
        network_data = WaterModels.parse_file(shamir_ts_tank_path)
        mn_data = WaterModels.make_multinetwork(network_data)
        wm = build_generic_model(mn_data, NCNLPWaterModel, WaterModels.post_mn_wf, multinetwork=true)
        solution = solve_generic_model(wm, ipopt)

        @test solution["termination_status"] == MOI.LOCALLY_SOLVED
        @test isapprox(solution["solution"]["nw"]["1"]["pipes"]["2"]["q"], 0.09356500, rtol=1.0e-4)
        @test isapprox(solution["solution"]["nw"]["1"]["pipes"]["6"]["q"], 0.05571000, rtol=1.0e-4)
        @test isapprox(solution["solution"]["nw"]["2"]["pipes"]["2"]["q"], 0.04678500, rtol=1.0e-4)
        @test isapprox(solution["solution"]["nw"]["2"]["pipes"]["6"]["q"], 0.02785300, rtol=1.0e-4)
        @test isapprox(solution["solution"]["nw"]["3"]["pipes"]["2"]["q"], 0.02339300, rtol=1.0e-4)
        @test isapprox(solution["solution"]["nw"]["3"]["pipes"]["6"]["q"], 0.01392600, rtol=1.0e-4)
        @test isapprox(solution["solution"]["nw"]["1"]["nodes"]["1"]["h"], 220.000000, rtol=1.0e-3)
        @test isapprox(solution["solution"]["nw"]["1"]["nodes"]["3"]["h"], 200.466599, rtol=1.0e-3)
        @test isapprox(solution["solution"]["nw"]["1"]["nodes"]["5"]["h"], 193.808380, rtol=1.0e-3)
        @test isapprox(solution["solution"]["nw"]["2"]["nodes"]["1"]["h"], 216.435196, rtol=1.0e-3)
        @test isapprox(solution["solution"]["nw"]["2"]["nodes"]["3"]["h"], 211.023941, rtol=1.0e-3)
        @test isapprox(solution["solution"]["nw"]["2"]["nodes"]["5"]["h"], 209.179306, rtol=1.0e-3)
        @test isapprox(solution["solution"]["nw"]["3"]["nodes"]["1"]["h"], 214.652771, rtol=1.0e-3)
        @test isapprox(solution["solution"]["nw"]["3"]["nodes"]["3"]["h"], 213.153824, rtol=1.0e-3)
        @test isapprox(solution["solution"]["nw"]["3"]["nodes"]["5"]["h"], 212.642853, rtol=1.0e-3)
    end

    @testset "Shamir network (unknown flow directions), MICP formulation." begin
        solution = run_wf(shamir_path, MICPWaterModel, ipopt, relaxed=true)
        @test solution["termination_status"] == MOI.LOCALLY_SOLVED
    end

    @testset "Shamir network (unknown flow directions), MILPR formulation." begin
        solution = run_wf(shamir_path, MILPRWaterModel, cbc)
        @test solution["termination_status"] == MOI.OPTIMAL
        @test isapprox(solution["solution"]["pipes"]["2"]["q"], 0.093565, rtol=0.05)
        @test isapprox(solution["solution"]["pipes"]["6"]["q"], 0.055710, rtol=0.05)
        @test isapprox(solution["solution"]["nodes"]["2"]["h"], 203.247650, rtol=0.05)
        @test isapprox(solution["solution"]["nodes"]["6"]["h"], 195.445953, rtol=0.05)
    end
end
