@testset "Hazen-Williams Non-convex MINLP Problems" begin
    #@testset "foss_poly_1" begin
    #    network_path = "../test/data/epanet/foss_poly_1.inp"
    #    solution = run_wf_hw(network_path, MINLPWaterModel, bonmin)
    #    @test solution["status"] == :LocalOptimal
    #end

    #@testset "hanoi" begin
    #    network_path = "../test/data/epanet/hanoi.inp"
    #    solution = run_wf_hw(network_path, MINLPWaterModel, bonmin)
    #    InfrastructureModels.print_summary(solution["solution"])
    #    @test solution["status"] == :LocalOptimal
    #end

    #@testset "hanoi_extended" begin
    #    network_path = "../test/data/epanet/hanoi_extended.inp"
    #    solution = run_wf_hw(network_path, MINLPWaterModel, bonmin)
    #    @test solution["status"] == :LocalOptimal
    #end

    #@testset "zj" begin
    #    network_path = "../test/data/epanet/zj.inp"
    #    solution = run_wf_hw(network_path, MINLPWaterModel, bonmin)
    #    @test solution["status"] == :LocalOptimal
    #end

    #@testset "klmod" begin
    #    network_path = "../test/data/epanet/klmod.inp"
    #    solution = run_wf_hw(network_path, MINLPWaterModel, bonmin)
    #    @test solution["status"] == :LocalOptimal
    #end
end

@testset "Hazen-Williams Non-convex MINLP Problems with Fixed Directions" begin
    #@testset "foss_poly_1" begin
    #    network_path = "../test/data/epanet/foss_poly_1.inp"
    #    solution = run_wf_hw(network_path, MINLPWaterModel, bonmin)
    #    @test solution["status"] == :LocalOptimal
    #end

    @testset "hanoi" begin
        network_path = "../test/data/epanet/hanoi.inp"
        modification_path = "../test/data/json/wf-hanoi.json"
        solution = run_wf_hw(network_path, modification_path, MINLPWaterModel, bonmin)
        InfrastructureModels.print_summary(solution["solution"])
        @test solution["status"] == :LocalOptimal
    end

    #@testset "hanoi_extended" begin
    #    network_path = "../test/data/epanet/hanoi_extended.inp"
    #    solution = run_wf_hw(network_path, MINLPWaterModel, bonmin)
    #    @test solution["status"] == :LocalOptimal
    #end

    #@testset "zj" begin
    #    network_path = "../test/data/epanet/zj.inp"
    #    solution = run_wf_hw(network_path, MINLPWaterModel, bonmin)
    #    @test solution["status"] == :LocalOptimal
    #end

    #@testset "klmod" begin
    #    network_path = "../test/data/epanet/klmod.inp"
    #    solution = run_wf_hw(network_path, MINLPWaterModel, bonmin)
    #    @test solution["status"] == :LocalOptimal
    #end
end

#@testset "Hazen-Williams Convex MINLP Problems" begin
#    @testset "foss_poly_1" begin
#        network_path = "../test/data/epanet/foss_poly_1.inp"
#        solution = run_wf_hw(network_path, ConvexMINLPWaterModel, pavito)
#        @test solution["status"] == :LocalOptimal
#    end
#
#    @testset "hanoi" begin
#        network_path = "../test/data/epanet/hanoi.inp"
#        solution = run_wf_hw(network_path, ConvexMINLPWaterModel, pavito)
#        @test solution["status"] == :LocalOptimal
#    end
#
#    @testset "hanoi_extended" begin
#        network_path = "../test/data/epanet/hanoi_extended.inp"
#        solution = run_wf_hw(network_path, ConvexMINLPWaterModel, pavito)
#        @test solution["status"] == :LocalOptimal
#    end
#
#    @testset "zj" begin
#        network_path = "../test/data/epanet/zj.inp"
#        solution = run_wf_hw(network_path, ConvexMINLPWaterModel, pavito)
#        @test solution["status"] == :LocalOptimal
#    end
#
#    @testset "klmod" begin
#        network_path = "../test/data/epanet/klmod.inp"
#        solution = run_wf_hw(network_path, ConvexMINLPWaterModel, pavito)
#        @test solution["status"] == :LocalOptimal
#    end
#end
#
#@testset "Hazen-Williams Standard MILP Problems" begin
#    #@testset "foss_poly_1" begin
#    #    network_path = "../test/data/epanet/foss_poly_1.inp"
#    #    solution = run_wf_hw(network_path, MILPWaterModel, cbc)
#    #    @test solution["status"] == :LocalOptimal
#    #end
#
#    @testset "hanoi" begin
#        network_path = "../test/data/epanet/hanoi.inp"
#        solution = run_wf_hw(network_path, MILPWaterModel, cbc)
#        @test solution["status"] == :Optimal
#    end
#
#    @testset "hanoi_extended" begin
#        network_path = "../test/data/epanet/hanoi_extended.inp"
#        solution = run_wf_hw(network_path, MILPWaterModel, cbc)
#        @test solution["status"] == :Optimal
#    end
#
#    #@testset "zj" begin
#    #    network_path = "../test/data/epanet/zj.inp"
#    #    solution = run_wf_hw(network_path, MILPWaterModel, cbc)
#    #    @test solution["status"] == :LocalOptimal
#    #end
#
#    #@testset "klmod" begin
#    #    network_path = "../test/data/epanet/klmod.inp"
#    #    solution = run_wf_hw(network_path, MILPWaterModel, cbc)
#    #    @test solution["status"] == :LocalOptimal
#    #end
#end
#
#@testset "Hazen-Williams Relaxed MILP Problems" begin
#    @testset "foss_poly_1" begin
#        network_path = "../test/data/epanet/foss_poly_1.inp"
#        solution = run_wf_hw(network_path, RelaxedMILPWaterModel, cbc)
#        @test solution["status"] == :Optimal
#    end
#
#    @testset "hanoi" begin
#        network_path = "../test/data/epanet/hanoi.inp"
#        solution = run_wf_hw(network_path, RelaxedMILPWaterModel, cbc)
#        @test solution["status"] == :Optimal
#    end
#
#    @testset "hanoi_extended" begin
#        network_path = "../test/data/epanet/hanoi_extended.inp"
#        solution = run_wf_hw(network_path, RelaxedMILPWaterModel, cbc)
#        @test solution["status"] == :Optimal
#    end
#
#    @testset "zj" begin
#        network_path = "../test/data/epanet/zj.inp"
#        solution = run_wf_hw(network_path, RelaxedMILPWaterModel, cbc)
#        @test solution["status"] == :Optimal
#    end
#
#    @testset "klmod" begin
#        network_path = "../test/data/epanet/klmod.inp"
#        solution = run_wf_hw(network_path, RelaxedMILPWaterModel, cbc)
#        @test solution["status"] == :Optimal
#    end
#end
#
#@testset "Darcy-Weisbach Non-convex MINLP Problems" begin
#    #@testset "balerma" begin
#    #    network_path = "../test/data/epanet/balerma.inp"
#    #    solution = run_wf_dw(network_path, MINLPWaterModel, bonmin)
#    #    @test solution["status"] == :LocalOptimal
#    #end
#
#    #@testset "rural" begin
#    #    network_path = "../test/data/epanet/rural.inp"
#    #    solution = run_wf_dw(network_path, MINLPWaterModel, bonmin)
#    #    @test solution["status"] == :LocalOptimal
#    #end
#end
#
#@testset "Darcy-Weisbach Convex MINLP Problems" begin
#    @testset "balerma" begin
#        network_path = "../test/data/epanet/balerma.inp"
#        solution = run_wf_dw(network_path, ConvexMINLPWaterModel, pavito)
#        @test solution["status"] == :LocalOptimal
#    end
#
#    @testset "rural" begin
#        network_path = "../test/data/epanet/rural.inp"
#        solution = run_wf_dw(network_path, ConvexMINLPWaterModel, pavito)
#        @test solution["status"] == :LocalOptimal
#    end
#end
#
#@testset "Darcy-Weisbach Standard MILP Problems" begin
#    #@testset "balerma" begin
#    #    network_path = "../test/data/epanet/balerma.inp"
#    #    solution = run_wf_dw(network_path, MILPWaterModel, cbc)
#    #    @test solution["status"] == :LocalOptimal
#    #end
#
#    #@testset "rural" begin
#    #    network_path = "../test/data/epanet/rural.inp"
#    #    solution = run_wf_dw(network_path, MILPWaterModel, cbc)
#    #    @test solution["status"] == :LocalOptimal
#    #end
#end
#
#@testset "Darcy-Weisbach Relaxed MILP Problems" begin
#    @testset "balerma" begin
#        network_path = "../test/data/epanet/balerma.inp"
#        solution = run_wf_dw(network_path, RelaxedMILPWaterModel, cbc)
#        @test solution["status"] == :Optimal
#    end
#
#    @testset "rural" begin
#        network_path = "../test/data/epanet/rural.inp"
#        solution = run_wf_dw(network_path, RelaxedMILPWaterModel, cbc)
#        @test solution["status"] == :Optimal
#    end
#end
