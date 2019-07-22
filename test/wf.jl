@testset "Water Flow Problems" begin
    # Set test network paths.
    balerma_path = "../test/data/epanet/balerma.inp"
    richmond_skeleton_path = "../test/data/epanet/richmond-skeleton.inp"
    shamir_path = "../test/data/epanet/shamir.inp"

    @testset "Balerma network (unknown flow directions), CNLP formulation." begin
        solution = run_wf(balerma_path, CNLPWaterModel, ipopt)
        @test solution["termination_status"] == MOI.LOCALLY_SOLVED
    end

    #@testset "Richmond network (unknown flow directions), NCNLP formulation." begin
    #    solution = run_wf(richmond_skeleton_path, NCNLPWaterModel, ipopt)
    #end

    @testset "Shamir network (unknown flow directions), CNLP formulation." begin
        shamir_mn_data = build_mn_data("../test/data/epanet/shamir.inp")
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
        @test isapprox(solution["solution"]["junctions"]["2"]["h"], 203.247650, rtol=1.0e-4)
        @test isapprox(solution["solution"]["junctions"]["6"]["h"], 195.445953, rtol=1.0e-4)
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
        @test isapprox(solution["solution"]["junctions"]["2"]["h"], 203.247650, rtol=0.05)
        @test isapprox(solution["solution"]["junctions"]["6"]["h"], 195.445953, rtol=0.05)
    end
end
