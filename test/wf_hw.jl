@testset "Hazen-Williams NLP Problems" begin
    @testset "Hanoi network (unknown flow directions)." begin
        network_path = "../test/data/epanet/hanoi.inp"
        solution = run_wf_hw(network_path, NLPWaterModel, ipopt)
        InfrastructureModels.print_summary(solution["solution"])
        @test solution["status"] == :LocalOptimal
    end

    @testset "Hanoi network (known flow directions)." begin
        network_path = "../test/data/epanet/hanoi.inp"
        modification_path = "../test/data/json/wf-hanoi.json"
        solution = run_wf_hw(network_path, modification_path, NLPWaterModel, ipopt)
        @test solution["status"] == :LocalOptimal
    end
end

#@testset "Hazen-Williams MINLP-B Problems" begin
#    @testset "Hanoi network (unknown flow directions)." begin
#        network_path = "../test/data/epanet/hanoi.inp"
#        solution = run_wf_hw(network_path, MINLPBWaterModel, bonmin)
#        @test solution["status"] == :LocalOptimal
#    end
#
#    @testset "Hanoi network (known flow directions)." begin
#        network_path = "../test/data/epanet/hanoi.inp"
#        modification_path = "../test/data/json/wf-hanoi.json"
#        solution = run_wf_hw(network_path, modification_path, MINLPBWaterModel, bonmin)
#        @test solution["status"] == :LocalOptimal
#    end
#end
#
#@testset "Hazen-Williams MICP Problems" begin
#    @testset "Hanoi network (unknown flow directions)." begin
#        network_path = "../test/data/epanet/hanoi.inp"
#        solution = run_wf_hw(network_path, MICPWaterModel, bonmin)
#        @test solution["status"] == :LocalOptimal
#    end
#
#    @testset "Hanoi network (known flow directions)." begin
#        network_path = "../test/data/epanet/hanoi.inp"
#        modification_path = "../test/data/json/wf-hanoi.json"
#        solution = run_wf_hw(network_path, modification_path, MICPWaterModel, bonmin)
#        @test solution["status"] == :LocalOptimal
#    end
#end
#
#@testset "Hazen-Williams MILP Problems" begin
#    @testset "Hanoi network (unknown flow directions)." begin
#        network_path = "../test/data/epanet/hanoi.inp"
#        solution = run_wf_hw(network_path, MILPWaterModel, cbc)
#        @test solution["status"] == :Optimal
#    end
#
#    @testset "Hanoi network (known flow directions)." begin
#        network_path = "../test/data/epanet/hanoi.inp"
#        modification_path = "../test/data/json/wf-hanoi.json"
#        solution = run_wf_hw(network_path, modification_path, MILPWaterModel, cbc)
#        @test solution["status"] == :Optimal
#    end
#end
#
#@testset "Hazen-Williams MILP-R Problems" begin
#    @testset "Hanoi network (unknown flow directions)." begin
#        network_path = "../test/data/epanet/hanoi.inp"
#        solution = run_wf_hw(network_path, MILPRWaterModel, cbc)
#        @test solution["status"] == :Optimal
#    end
#
#    @testset "Hanoi network (known flow directions)." begin
#        network_path = "../test/data/epanet/hanoi.inp"
#        modification_path = "../test/data/json/wf-hanoi.json"
#        solution = run_wf_hw(network_path, modification_path, MILPRWaterModel, cbc)
#        @test solution["status"] == :Optimal
#    end
#end
