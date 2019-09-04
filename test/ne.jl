@testset "Network Expansion Problems" begin
    # Get an initial solution for the Hazen-Williams network.
    hw_network_path = "../test/data/epanet/shamir.inp"
    hw_network = WaterModels.parse_file(hw_network_path)
    hw_wf_result = WaterModels.run_wf(hw_network, NCNLPWaterModel, ipopt)
    InfrastructureModels.update_data!(hw_network, hw_wf_result["solution"])

    #@testset "Shamir network (unknown flow directions), CNLP formulation." begin
    #    network = deepcopy(hw_network)
    #    modifications = WaterModels.parse_file("../test/data/json/shamir.json")
    #    InfrastructureModels.update_data!(network, modifications)
    #    @test_throws ErrorException run_ne(network, CNLPWaterModel, ipopt, relaxed=true)
    #end

    @testset "Shamir network (unknown flow directions), MICP formulation." begin
        network = WaterModels.parse_file(hw_network_path)
        modifications = WaterModels.parse_file("../test/data/json/shamir-reduced.json")
        InfrastructureModels.update_data!(network, modifications)
        wm = build_model(network, MICPWaterModel, WaterModels.post_ne, multinetwork=false)
        f = Juniper.register(fun(wm, :head_loss)..., autodiff=false)
        juniper = JuMP.with_optimizer(Juniper.Optimizer, nl_solver=ipopt, registered_functions=[f], log_levels=[])
        solution = WaterModels.optimize_model!(wm, juniper)
        @test solution["termination_status"] == MOI.LOCALLY_SOLVED
    end

    @testset "Shamir network (unknown flow directions), MILPR formulation." begin
        network = WaterModels.parse_file(hw_network_path)
        modifications = WaterModels.parse_file("../test/data/json/shamir-reduced.json")
        InfrastructureModels.update_data!(network, modifications)
        result = run_ne(network, MILPRWaterModel, cbc)
        @test result["objective_value"] == 1.36e6
        @test result["termination_status"] == MOI.OPTIMAL
    end

    @testset "Shamir network (unknown flow directions), NCNLP formulation." begin
        network = WaterModels.parse_file(hw_network_path)
        modifications = WaterModels.parse_file("../test/data/json/shamir-reduced.json")
        InfrastructureModels.update_data!(network, modifications)
        wm = build_model(network, NCNLPWaterModel, WaterModels.post_ne, multinetwork=false)
        f = Juniper.register(fun(wm, :head_loss)..., autodiff=false)
        juniper = JuMP.with_optimizer(Juniper.Optimizer, nl_solver=ipopt, registered_functions=[f], log_levels=[])
        solution = WaterModels.optimize_model!(wm, juniper)
        @test solution["termination_status"] == MOI.LOCALLY_SOLVED
    end
end
