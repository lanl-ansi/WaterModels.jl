@testset "Darcy-Weisbach NLP Problems" begin
    @testset "Balerma network (unknown flow directions)." begin
        network_path = "../test/data/epanet/balerma.inp"
        solution = run_wf_dw(network_path, NLPWaterModel, ipopt)
        @test solution["status"] == :LocalOptimal
    end

    @testset "Balerma network (known flow directions)." begin
        network_path = "../test/data/epanet/balerma.inp"
        modification_path = "../test/data/json/wf-balerma.json"
        solution = run_wf_dw(network_path, modification_path, NLPWaterModel, ipopt)
        @test solution["status"] == :LocalOptimal
    end
end

@testset "Darcy-Weisbach MICP Problems" begin
    #@testset "Balerma network (unknown flow directions)." begin
    #    network_path = "../test/data/epanet/balerma.inp"
    #    solution = run_wf_dw(network_path, MICPWaterModel, bonmin)
    #    @test solution["status"] == :LocalOptimal
    #end

    @testset "Balerma network (known flow directions)." begin
        network_path = "../test/data/epanet/balerma.inp"
        modification_path = "../test/data/json/wf-balerma.json"
        solution = run_wf_dw(network_path, modification_path, MICPWaterModel, bonmin)
        @test solution["status"] == :LocalOptimal
    end
end

@testset "Darcy-Weisbach MILP Problems" begin
    #@testset "Balerma network (unknown flow directions)." begin
    #    network_path = "../test/data/epanet/balerma.inp"
    #    solution = run_wf_dw(network_path, MILPWaterModel, cbc)
    #    @test solution["status"] == :Optimal
    #end

    #@testset "Balerma network (known flow directions)." begin
    #    network_path = "../test/data/epanet/balerma.inp"
    #    modification_path = "../test/data/json/wf-balerma.json"
    #    solution = run_wf_dw(network_path, modification_path, MILPWaterModel, cbc)
    #    @test solution["status"] == :Optimal
    #end
end

@testset "Darcy-Weisbach MILP-R Problems" begin
    @testset "Balerma network (unknown flow directions)." begin
        network_path = "../test/data/epanet/balerma.inp"
        solution = run_wf_dw(network_path, MILPRWaterModel, cbc)
        @test solution["status"] == :Optimal
    end

    @testset "Balerma network (known flow directions)." begin
        network_path = "../test/data/epanet/balerma.inp"
        modification_path = "../test/data/json/wf-balerma.json"
        solution = run_wf_dw(network_path, modification_path, MILPRWaterModel, cbc)
        @test solution["status"] == :Optimal
    end
end
