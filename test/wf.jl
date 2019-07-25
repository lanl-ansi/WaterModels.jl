@testset "Water Flow Problems" begin
    # Set test network paths.
    balerma_path = "../test/data/epanet/balerma.inp"
    richmond_skeleton_path = "../test/data/epanet/richmond-skeleton.inp"
    richmond_skeleton_sp_path = "../test/data/epanet/richmond-skeleton-sp.inp"
    shamir_path = "../test/data/epanet/shamir.inp"
    shamir_ts_path = "../test/data/epanet/shamir-ts.inp"

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
        @test isapprox(solution["solution"]["nw"]["1"]["pipes"]["2"]["q"], 0.093565, rtol=1.0e-4)
        @test isapprox(solution["solution"]["nw"]["1"]["pipes"]["6"]["q"], 0.055710, rtol=1.0e-4)
        @test isapprox(solution["solution"]["nw"]["2"]["pipes"]["2"]["q"], 0.046785, rtol=1.0e-4)
        @test isapprox(solution["solution"]["nw"]["2"]["pipes"]["6"]["q"], 0.027853, rtol=1.0e-4)
        @test isapprox(solution["solution"]["nw"]["3"]["pipes"]["2"]["q"], 0.023393, rtol=1.0e-4)
        @test isapprox(solution["solution"]["nw"]["3"]["pipes"]["6"]["q"], 0.013926, rtol=1.0e-4)
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
