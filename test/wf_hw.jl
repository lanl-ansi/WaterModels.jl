@testset "Hazen-Williams NLP Problems" begin
    @testset "Shamir network (unknown flow directions)." begin
        network_path = "../test/data/epanet/shamir.inp"
        solution = run_wf_hw(network_path, NLPWaterModel, ipopt)
        @test solution["status"] == :LocalOptimal
    end

    #@testset "Shamir network (known flow directions)." begin
    #    network_path = "../test/data/epanet/shamir.inp"
    #    modification_path = "../test/data/json/shamir.json"
    #    solution = run_wf_hw(network_path, modification_path, NLPWaterModel, ipopt)
    #    @test solution["status"] == :LocalOptimal
    #end
end

@testset "Hazen-Williams MINLP-B Problems" begin
    #@testset "Shamir network (unknown flow directions)." begin
    #    network_path = "../test/data/epanet/shamir.inp"
    #    solution = run_wf_hw(network_path, MINLPBWaterModel, bonmin)
    #    @test solution["status"] == :LocalOptimal
    #end

    #@testset "Shamir network (known flow directions)." begin
    #    network_path = "../test/data/epanet/shamir.inp"
    #    modification_path = "../test/data/json/shamir.json"
    #    solution = run_wf_hw(network_path, modification_path, MINLPBWaterModel, bonmin)
    #    @test solution["status"] == :LocalOptimal
    #end
end

@testset "Hazen-Williams MICP Problems" begin
    #@testset "Shamir network (unknown flow directions)." begin
    #    network_path = "../test/data/epanet/shamir.inp"
    #    solution = run_wf_hw(network_path, MICPWaterModel, bonmin)
    #    @test solution["status"] == :LocalOptimal
    #end

    #@testset "Shamir network (known flow directions)." begin
    #    network_path = "../test/data/epanet/shamir.inp"
    #    modification_path = "../test/data/json/shamir.json"
    #    solution = run_wf_hw(network_path, modification_path, MICPWaterModel, bonmin)
    #    @test solution["status"] == :LocalOptimal
    #end
end

@testset "Hazen-Williams MILP Problems" begin
    #@testset "Shamir network (unknown flow directions)." begin
    #    network_path = "../test/data/epanet/shamir.inp"
    #    solution = run_wf_hw(network_path, MILPWaterModel, cbc)
    #    @test solution["status"] == :Optimal
    #end

    #@testset "Shamir network (known flow directions)." begin
    #    network_path = "../test/data/epanet/shamir.inp"
    #    modification_path = "../test/data/json/shamir.json"
    #    solution = run_wf_hw(network_path, modification_path, MILPWaterModel, cbc)
    #    @test solution["status"] == :Optimal
    #end
end

@testset "Hazen-Williams MILP-R Problems" begin
    @testset "Shamir network (unknown flow directions)." begin
        network_path = "../test/data/epanet/shamir.inp"
        solution = run_wf_hw(network_path, MILPRWaterModel, cbc)
        @test solution["status"] == :Optimal
    end

    #@testset "Shamir network (known flow directions)." begin
    #    network_path = "../test/data/epanet/shamir.inp"
    #    modification_path = "../test/data/json/shamir.json"
    #    solution = run_wf_hw(network_path, modification_path, MILPRWaterModel, cbc)
    #    @test solution["status"] == :Optimal
    #end
end
