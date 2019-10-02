@testset "Network Expansion Problems" begin
    @testset "Shamir network (reduced), CNLP formulation." begin
        data = WaterModels.parse_file("../test/data/epanet/shamir.inp")
        modifications = WaterModels.parse_file("../test/data/json/shamir-reduced.json")
        InfrastructureModels.update_data!(data, modifications)
        @test_throws ErrorException build_model(data, CNLPWaterModel, WaterModels.post_ne)
    end

    @testset "Shamir network (reduced), MICP formulation." begin
        data = WaterModels.parse_file("../test/data/epanet/shamir.inp")
        modifications = WaterModels.parse_file("../test/data/json/shamir-reduced.json")
        InfrastructureModels.update_data!(data, modifications)
        wm = build_model(data, MICPWaterModel, WaterModels.post_ne)
        f = Juniper.register(fun(wm, :head_loss)..., autodiff=false)
        juniper = JuMP.with_optimizer(Juniper.Optimizer, nl_solver=ipopt, registered_functions=[f], log_levels=[])
        solution = WaterModels.optimize_model!(wm, juniper)

        @test solution["termination_status"] == LOCALLY_SOLVED
        @test isapprox(solution["objective"], 1.36e6, rtol=1.0e-4)
    end

    #@testset "Shamir network (unknown flow directions), MILPR formulation." begin
    #    network = WaterModels.parse_file(hw_network_path)
    #    modifications = WaterModels.parse_file("../test/data/json/shamir-reduced.json")
    #    InfrastructureModels.update_data!(network, modifications)
    #    result = run_ne(network, MILPRWaterModel, cbc)
    #    @test result["objective"] == 1.36e6
    #    @test result["termination_status"] == OPTIMAL
    #end
    
    @testset "Shamir network (reduced), NCNLP formulation." begin
        data = WaterModels.parse_file("../test/data/epanet/shamir.inp")
        modifications = WaterModels.parse_file("../test/data/json/shamir-reduced.json")
        InfrastructureModels.update_data!(data, modifications)
        wm = build_model(data, NCNLPWaterModel, WaterModels.post_ne)
        f = Juniper.register(fun(wm, :head_loss)..., autodiff=false)
        juniper = JuMP.with_optimizer(Juniper.Optimizer, nl_solver=ipopt, registered_functions=[f], log_levels=[])
        solution = WaterModels.optimize_model!(wm, juniper)

        @test solution["termination_status"] == LOCALLY_SOLVED
        @test isapprox(solution["objective"], 1.36e6, rtol=1.0e-4)
    end
end
