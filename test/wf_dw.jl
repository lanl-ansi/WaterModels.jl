@testset "Darcy-Weisbach Non-convex MINLP Problems" begin
    @testset "balerma" begin
        network_path = "../test/data/epanet/balerma.inp"
        solution = run_wf_dw(network_path, MINLPWaterModel, bonmin)
        @test solution["status"] == :LocalOptimal
    end
end

@testset "Darcy-Weisbach Non-convex MINLP Problems with Fixed Directions" begin
    @testset "balerma" begin
        network_path = "../test/data/epanet/balerma.inp"
        modification_path = "../test/data/json/wf-balerma.json"
        solution = run_wf_dw(network_path, modification_path, MINLPWaterModel, bonmin)
        @test solution["status"] == :LocalOptimal
    end
end

#@testset "Darcy-Weisbach MICP Problems" begin
#    @testset "balerma" begin
#        network_path = "../test/data/epanet/balerma.inp"
#        solution = run_wf_dw(network_path, MICPWaterModel, pavito)
#        @test solution["status"] == :Optimal
#    end
#end

@testset "Darcy-Weisbach MICP Problems with Fixed Directions" begin
    @testset "balerma" begin
        network_path = "../test/data/epanet/balerma.inp"
        modification_path = "../test/data/json/wf-balerma.json"
        solution = run_wf_dw(network_path, modification_path, MICPWaterModel, pavito)
        @test solution["status"] == :Optimal
    end
end

#@testset "Darcy-Weisbach Standard MILP Problems" begin
#    @testset "balerma" begin
#        network_path = "../test/data/epanet/balerma.inp"
#        solution = run_wf_dw(network_path, MILPWaterModel, cbc)
#        @test solution["status"] == :Optimal
#    end
#end

#@testset "Darcy-Weisbach Standard MILP Problems with Fixed Directions" begin
#    @testset "balerma" begin
#        network_path = "../test/data/epanet/balerma.inp"
#        modification_path = "../test/data/json/wf-balerma.json"
#        solution = run_wf_dw(network_path, modification_path, MILPWaterModel, cbc)
#        @test solution["status"] == :Optimal
#    end
#end

#@testset "Darcy-Weisbach Relaxed MILP Problems" begin
#    @testset "balerma" begin
#        network_path = "../test/data/epanet/balerma.inp"
#        solution = run_wf_dw(network_path, RelaxedMILPWaterModel, cbc)
#        @test solution["status"] == :Optimal
#    end
#end

@testset "Darcy-Weisbach Relaxed MILP Problems with Fixed Directions" begin
    @testset "balerma" begin
        network_path = "../test/data/epanet/balerma.inp"
        modification_path = "../test/data/json/wf-balerma.json"
        solution = run_wf_dw(network_path, modification_path, RelaxedMILPWaterModel, cbc)
        @test solution["status"] == :Optimal
    end
end
