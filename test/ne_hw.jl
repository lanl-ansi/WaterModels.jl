#@testset "Hazen-Williams NLP Problems" begin
#    @testset "Hanoi network." begin
#        network_path = "../test/data/epanet/hanoi.inp"
#        modification_path = "../test/data/json/ne-hanoi.json"
#        solution = run_ne_hw(network_path, modification_path, NLPWaterModel, bonmin)
#        InfrastructureModels.print_summary(solution["solution"])
#        @test solution["status"] == :LocalOptimal
#    end
#end

#@testset "Hazen-Williams MINLP-B Problems" begin
#    @testset "Hanoi network." begin
#        network_path = "../test/data/epanet/hanoi.inp"
#        modification_path = "../test/data/json/ne-hanoi.json"
#        solution = run_ne_hw(network_path, modification_path, MINLPBWaterModel, bonmin)
#        InfrastructureModels.print_summary(solution["solution"])
#        @test solution["status"] == :LocalOptimal
#    end
#end

@testset "Hazen-Williams MICP Problems" begin
    @testset "Hanoi network." begin
        network_path = "../test/data/epanet/hanoi.inp"
        modification_path = "../test/data/json/ne-hanoi.json"
        solution = run_ne_hw(network_path, modification_path, MICPWaterModel, bonmin)
        InfrastructureModels.print_summary(solution["solution"])
        @test solution["status"] == :LocalOptimal
    end
end

#@testset "Hazen-Williams MILP Problems" begin
#    @testset "Hanoi network." begin
#        network_path = "../test/data/epanet/hanoi.inp"
#        modification_path = "../test/data/json/ne-hanoi.json"
#        solution = run_ne_hw(network_path, modification_path, MILPWaterModel, cbc)
#        InfrastructureModels.print_summary(solution["solution"])
#        @test solution["status"] == :Optimal
#    end
#end

#@testset "Hazen-Williams MILP-R Problems" begin
#    @testset "Hanoi network." begin
#        network_path = "../test/data/epanet/hanoi.inp"
#        modification_path = "../test/data/json/ne-hanoi.json"
#        solution = run_ne_hw(network_path, modification_path, MILPRWaterModel, cbc)
#        InfrastructureModels.print_summary(solution["solution"])
#        @test solution["status"] == :Optimal
#    end
#end
