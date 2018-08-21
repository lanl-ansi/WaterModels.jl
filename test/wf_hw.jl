@testset "Hazen-Williams Non-convex MINLP Problems" begin
    @testset "hanoi" begin
        network_path = "../test/data/epanet/hanoi.inp"
        solution = run_wf_hw(network_path, MINLPWaterModel, bonmin)
        InfrastructureModels.print_summary(solution["solution"])
        @test solution["status"] == :LocalOptimal
    end
end

@testset "Hazen-Williams Non-convex MINLP Problems with Fixed Directions" begin
    @testset "hanoi" begin
        network_path = "../test/data/epanet/hanoi.inp"
        modification_path = "../test/data/json/wf-hanoi.json"
        solution = run_wf_hw(network_path, modification_path, MINLPWaterModel, bonmin)
        InfrastructureModels.print_summary(solution["solution"])
        @test solution["status"] == :LocalOptimal
    end
end

@testset "Hazen-Williams MICP Problems" begin
    @testset "hanoi" begin
        network_path = "../test/data/epanet/hanoi.inp"
        solution = run_wf_hw(network_path, MICPWaterModel, pavito)
        @test solution["status"] == :Optimal
    end
end

@testset "Hazen-Williams MICP Problems with Fixed Directions" begin
    @testset "hanoi" begin
        network_path = "../test/data/epanet/hanoi.inp"
        modification_path = "../test/data/json/wf-hanoi.json"
        solution = run_wf_hw(network_path, modification_path, MICPWaterModel, pavito)
        @test solution["status"] == :Optimal
    end
end

@testset "Hazen-Williams Standard MILP Problems" begin
    @testset "hanoi" begin
        network_path = "../test/data/epanet/hanoi.inp"
        solution = run_wf_hw(network_path, MILPWaterModel, cbc)
        @test solution["status"] == :Optimal
    end
end

@testset "Hazen-Williams Standard MILP Problems with Fixed Directions" begin
    @testset "hanoi" begin
        network_path = "../test/data/epanet/hanoi.inp"
        modification_path = "../test/data/json/wf-hanoi.json"
        solution = run_wf_hw(network_path, modification_path, MILPWaterModel, cbc)
        @test solution["status"] == :Optimal
    end
end

@testset "Hazen-Williams Relaxed MILP Problems" begin
    @testset "hanoi" begin
        network_path = "../test/data/epanet/hanoi.inp"
        solution = run_wf_hw(network_path, RelaxedMILPWaterModel, cbc)
        @test solution["status"] == :Optimal
    end
end

@testset "Hazen-Williams Relaxed MILP Problems with Fixed Directions" begin
    @testset "hanoi" begin
        network_path = "../test/data/epanet/hanoi.inp"
        modification_path = "../test/data/json/wf-hanoi.json"
        solution = run_wf_hw(network_path, modification_path, RelaxedMILPWaterModel, cbc)
        @test solution["status"] == :Optimal
    end
end
