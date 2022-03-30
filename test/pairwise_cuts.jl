@testset "src/util/pairwise_cuts.jl" begin
    data = WaterModels.parse_file("../test/data/epanet/multinetwork/owf-hw-lps.inp")
    mn_data = WaterModels.make_multinetwork(data)
    set_flow_partitions_si!(mn_data, 10.0, 1.0e-4)

    @testset "_PairwiseProblem instantiation" begin
        wm = instantiate_model(mn_data, LRDWaterModel, build_mn_wf)
        JuMP.set_optimizer(wm.model, milp_solver) # Explicitly set an optimizer.
        vid_1 = WaterModels._VariableIndex(1, :pump, :z_pump, 1)
        vid_2 = WaterModels._VariableIndex(1, :pipe, :y_pipe, 2)
        problem = WaterModels._PairwiseProblem(JuMP.MIN_SENSE, vid_1, vid_2, 0.0, 1.0e-4)

        @test problem.sense === JuMP.MIN_SENSE
        @test problem.variable_index_1 == vid_1
        @test problem.variable_index_2 == vid_2
        @test problem.variable_2_fixing_value == 0.0
    end

    @testset "_optimize_bound_problem!" begin
        wm = instantiate_model(mn_data, LRDWaterModel, build_mn_wf)
        JuMP.set_optimizer(wm.model, milp_solver) # Explicitly set an optimizer.
        vid_1 = WaterModels._VariableIndex(1, :pump, :z_pump, 1)
        vid_2 = WaterModels._VariableIndex(1, :pipe, :y_pipe, 2)
        problem = WaterModels._PairwiseProblem(JuMP.MIN_SENSE, vid_1, vid_2, 0.0, 1.0e-4)
        termination_status = WaterModels._optimize_bound_problem!(wm, problem)
        @test termination_status === OPTIMAL
    end

    @testset "_get_bound_problem_candidate! (minimization)" begin
        wm = instantiate_model(mn_data, LRDWaterModel, build_mn_wf)
        JuMP.set_optimizer(wm.model, milp_solver) # Explicitly set an optimizer.
        vid_1 = WaterModels._VariableIndex(1, :pump, :z_pump, 1)
        vid_2 = WaterModels._VariableIndex(1, :pipe, :y_pipe, 2)
        problem = WaterModels._PairwiseProblem(JuMP.MIN_SENSE, vid_1, vid_2, 0.0, 1.0e-4)
        @test WaterModels._get_bound_problem_candidate(wm, problem) == 0.0
    end

    @testset "_get_bound_problem_candidate! (maximization)" begin
        wm = instantiate_model(mn_data, LRDWaterModel, build_mn_wf)
        JuMP.set_optimizer(wm.model, milp_solver) # Explicitly set an optimizer.
        vid_1 = WaterModels._VariableIndex(1, :pump, :z_pump, 1)
        vid_2 = WaterModels._VariableIndex(1, :pipe, :y_pump, 1)
        problem = WaterModels._PairwiseProblem(JuMP.MAX_SENSE, vid_1, vid_2, 1.0, 1.0e-4)
        @test WaterModels._get_bound_problem_candidate(wm, problem) == 1.0
    end

    @testset "_get_bound_problem_candidate! (maximization, no solution)" begin
        wm = instantiate_model(mn_data, LRDWaterModel, build_mn_wf)
        JuMP.set_optimizer(wm.model, milp_solver) # Explicitly set an optimizer.
        vid_1 = WaterModels._VariableIndex(1, :pipe, :y_pipe, 2)
        vid_2 = WaterModels._VariableIndex(1, :pump, :z_pump, 1)
        problem = WaterModels._PairwiseProblem(JuMP.MAX_SENSE, vid_1, vid_2, 1.0, 1.0e-4)
        @test WaterModels._get_bound_problem_candidate(wm, problem) == 1.0
    end

    @testset "_unfix_bound_problem_variable!" begin
        wm = instantiate_model(mn_data, LRDWaterModel, build_mn_wf)
        JuMP.set_optimizer(wm.model, milp_solver) # Explicitly set an optimizer.
        vid_1 = WaterModels._VariableIndex(1, :pump, :z_pump, 1)
        vid_2 = WaterModels._VariableIndex(1, :pipe, :y_pipe, 2)
        problem = WaterModels._PairwiseProblem(JuMP.MIN_SENSE, vid_1, vid_2, 0.0, 1.0e-4)
        termination_status = WaterModels._optimize_bound_problem!(wm, problem)
        WaterModels._unfix_bound_problem_variable!(wm, problem)
        @test JuMP.is_fixed(WaterModels.var(wm, 1, :y_pipe, 2)) == false
    end

    @testset "_solve_bound_problem!" begin
        wm = instantiate_model(mn_data, LRDWaterModel, build_mn_wf)
        JuMP.set_optimizer(wm.model, milp_solver) # Explicitly set an optimizer.
        vid_1 = WaterModels._VariableIndex(1, :pump, :z_pump, 1)
        vid_2 = WaterModels._VariableIndex(1, :pipe, :y_pipe, 2)
        problem = WaterModels._PairwiseProblem(JuMP.MIN_SENSE, vid_1, vid_2, 0.0, 1.0e-4)
        @test WaterModels._solve_bound_problem!(wm, problem) == 0.0
    end

    @testset "_add_pairwise_cuts!(wm; nw = 1)" begin
        wm = instantiate_model(mn_data, LRDWaterModel, build_mn_wf)
        JuMP.set_optimizer(wm.model, milp_solver) # Explicitly set an optimizer.
        WaterModels._add_pairwise_cuts!(wm; nw = 1)
        result = WaterModels.optimize_model!(wm, optimizer=milp_solver)
        @test result["termination_status"] == OPTIMAL
    end
end
