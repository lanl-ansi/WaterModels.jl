@testset "Network Expansion Problems" begin
    # Get an initial solution for the Hazen-Williams network.
    hw_network_path = "../test/data/epanet/shamir.inp"
    hw_network = WaterModels.parse_file(hw_network_path)
    hw_wf_result = WaterModels.run_wf(hw_network, NCNLPWaterModel, ipopt, alpha=1.852)
    InfrastructureModels.update_data!(hw_network, hw_wf_result["solution"])

    @testset "Shamir network (unknown flow directions), CNLP formulation." begin
        network = deepcopy(hw_network)
        modifications = WaterModels.parse_file("../test/data/json/shamir.json")
        InfrastructureModels.update_data!(network, modifications)
        @test_throws ErrorException run_ne(network, CNLPWaterModel, ipopt, alpha=1.852, relaxed=true)
    end

    @testset "Shamir network (unknown flow directions), MICP formulation." begin
        network = deepcopy(hw_network)
        modifications = WaterModels.parse_file("../test/data/json/shamir-reduced.json")
        InfrastructureModels.update_data!(network, modifications)

        ## TODO: Use Juniper once user-defined derivatives are allowed.
        #post_ne = WaterModels.get_post_ne(1.852)
        #wm = build_generic_model(network, MICPWaterModel, post_ne)
        #f_1, f_2, f_3, f_4, f_5 = WaterModels.function_f_alpha_args(wm)
        #f = Juniper.register(f_1, f_2, f_3, f_4, f_5, autodiff=false)
        #juniper = JuMP.with_optimizer(Juniper.Optimizer, nl_solver=ipopt, registered_functions=[f])
        #result = WaterModels.solve_generic_model(wm, juniper)

        result = run_ne(network, MICPWaterModel, ipopt, alpha=1.852, relaxed=true)
        @test result["termination_status"] == MOI.ALMOST_LOCALLY_SOLVED ||
              result["termination_status"] == MOI.LOCALLY_SOLVED
    end

    @testset "Shamir network (unknown flow directions), MILPR formulation." begin
        network = deepcopy(hw_network)
        modifications = WaterModels.parse_file("../test/data/json/shamir-reduced.json")
        InfrastructureModels.update_data!(network, modifications)
        result = run_ne(network, MILPRWaterModel, cbc, alpha=1.852)
        @test result["objective_value"] == 1.36e6
        @test result["termination_status"] == MOI.OPTIMAL
    end

    @testset "Shamir network (unknown flow directions), NCNLP formulation." begin
        network = deepcopy(hw_network)
        modifications = WaterModels.parse_file("../test/data/json/shamir.json")
        InfrastructureModels.update_data!(network, modifications)
        result = run_ne(network, NCNLPWaterModel, ipopt_ws, alpha=1.852, relaxed=true)
        @test result["termination_status"] == MOI.ALMOST_LOCALLY_SOLVED ||
              result["termination_status"] == MOI.LOCALLY_SOLVED
    end
end
