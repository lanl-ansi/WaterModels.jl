@testset "Hazen-Williams MINLP Problems" begin
    #@testset "Shamir network (physical feasibility, relaxation)." begin
    #    network_path = "../test/data/epanet/shamir.inp"
    #    solution = run_wf_hw(network_path, MINLPWaterModel, pavito)
    #    status = solution["status"]
    #    @test status == :LocalOptimal || status == :Optimal
    #end

    #@testset "Shamir network (diameter selection (reduced), relaxation)." begin
    #    network_path = "../test/data/epanet/shamir.inp"
    #    modification_path = "../test/data/json/shamir-reduced.json"
    #    solution = run_ne_hw(network_path, modification_path, MINLPWaterModel, pavito)
    #    status = solution["status"]
    #    @test status == :LocalOptimal || status == :Optimal
    #end

    #@testset "Shamir network (diameter selection (reduced), global algorithm)." begin
    #    network_path = "../test/data/epanet/shamir.inp"
    #    modification_path = "../test/data/json/shamir-reduced.json"
    #    status = solve_global(network_path, modification_path, ipopt, glpk)
    #    @test status == :LocalOptimal || status == :Optimal
    #end

    @testset "Shamir network (diameter selection, global algorithm)." begin
        network_path = "../test/data/epanet/shamir.inp"
        modification_path = "../test/data/json/shamir.json"
        status = solve_global(network_path, modification_path, ipopt, gurobi)
        @test status == :LocalOptimal || status == :Optimal
    end
end

#@testset "Hazen-Williams NLP Problems" begin
#    @testset "Hanoi network." begin
#        network_path = "../test/data/epanet/shamir.inp"
#        modification_path = "../test/data/json/ne-shamir.json"
#        solution = run_ne_hw(network_path, modification_path, NLPWaterModel, bonmin)
#        InfrastructureModels.print_summary(solution["solution"])
#        @test solution["status"] == :LocalOptimal
#    end
#end
#
#@testset "Hazen-Williams MINLP-B Problems" begin
#    @testset "Hanoi network." begin
#        network_path = "../test/data/epanet/shamir.inp"
#        modification_path = "../test/data/json/ne-shamir.json"
#        solution = run_ne_hw(network_path, modification_path, MINLPBWaterModel, bonmin)
#        InfrastructureModels.print_summary(solution["solution"])
#        @test solution["status"] == :LocalOptimal
#    end
#end

#@testset "Hazen-Williams MICP Problems" begin
#    @testset "Hanoi network." begin
#        network_path = "../test/data/epanet/shamir.inp"
#        modification_path = "../test/data/json/ne-shamir.json"
#        solution = run_ne_hw(network_path, modification_path, MICPWaterModel, bonmin)
#        InfrastructureModels.print_summary(solution["solution"])
#        @test solution["status"] == :LocalOptimal
#    end
#end

#@testset "Hazen-Williams MILP Problems" begin
#    @testset "Hanoi network." begin
#        network_path = "../test/data/epanet/shamir.inp"
#        modification_path = "../test/data/json/ne-shamir.json"
#        solution = run_ne_hw(network_path, modification_path, MILPWaterModel, cbc)
#        InfrastructureModels.print_summary(solution["solution"])
#        @test solution["status"] == :Optimal
#    end
#end

#@testset "Hazen-Williams MILP-R Problems" begin
#    @testset "Hanoi network." begin
#        network_path = "../test/data/epanet/shamir.inp"
#        modification_path = "../test/data/json/shamir.json"
#        solution = run_ne_hw(network_path, modification_path, MILPRWaterModel, cbc)
#        @test solution["status"] == :Optimal
#    end
#end
