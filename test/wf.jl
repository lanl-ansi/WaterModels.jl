function get_hw_solution(network_path::String)
    model_type = GenericWaterModel{ConvexMINLPForm}
    solution = run_wf_hw(network_path, model_type, bonmin)
    InfrastructureModels.print_summary(solution["solution"])
    return solution
end

function get_dw_solution(network_path::String)
    model_type = GenericWaterModel{StandardMINLPForm}
    solution = run_wf_dw(network_path, model_type, bonmin)
    InfrastructureModels.print_summary(solution["solution"])
    return solution
end

@testset "Hazen-Williams MINLP Problems" begin
    #@testset "foss_poly_1" begin
    #    solution = get_hw_solution("../test/data/epanet/foss_poly_1.inp")
    #    @test solution["status"] == :LocalOptimal
    #end

    #@testset "hanoi_extended" begin
    #    solution = get_hw_solution("../test/data/epanet/hanoi_extended.inp")
    #    @test solution["status"] == :LocalOptimal
    #end

    @testset "hanoi" begin
        solution = get_hw_solution("../test/data/epanet/hanoi.inp")
        @test solution["status"] == :LocalOptimal
    end

    #@testset "zj" begin
    #    solution = get_hw_solution("../test/data/epanet/zj.inp")
    #    @test solution["status"] == :LocalOptimal
    #end

    # This one is too large to solve in a unit test using Bonmin.
    # @testset "klmod" begin
    #     solution = get_hw_solution("../test/data/epanet/klmod.inp")
    #     @test solution["status"] == :LocalOptimal
    # end
end

@testset "Darcy-Weisbach MINLP Problems" begin
    #@testset "balerma" begin
    #    solution = get_dw_solution("../test/data/epanet/balerma.inp")
    #    @test solution["status"] == :LocalOptimal
    #end

    # # TODO: Why doesn't a feasible solution exist to this problem?
    # @testset "rural" begin
    #     solution = get_dw_solution("../test/data/epanet/rural.inp")
    #     @test solution["status"] == :LocalOptimal
    # end
end
