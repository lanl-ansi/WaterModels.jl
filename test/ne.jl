@testset "Network Expansion Problems" begin
    #@testset "Shamir network (reduced), CNLP formulation." begin
    #    data = WaterModels.parse_file("../test/data/epanet/shamir.inp")
    #    modifications = WaterModels.parse_file("../test/data/json/shamir-reduced.json")
    #    InfrastructureModels.update_data!(data, modifications)
    #    @test_throws ErrorException build_model(data, CNLPWaterModel, WaterModels.post_ne)
    #end

    #@testset "Shamir network (reduced), MICP formulation." begin
    #    data = WaterModels.parse_file("../test/data/epanet/shamir.inp")
    #    modifications = WaterModels.parse_file("../test/data/json/shamir.json")
    #    InfrastructureModels.update_data!(data, modifications)
    #    wm = build_model(data, MICPWaterModel, WaterModels.post_ne)
    #    f = Juniper.register(fun(wm, :head_loss)..., autodiff=false)
    #    juniper = JuMP.with_optimizer(Juniper.Optimizer, nl_solver=ipopt, registered_functions=[f]) #, log_levels=[])
    #    solution = WaterModels.optimize_model!(wm, juniper)

    #    @test solution["termination_status"] == LOCALLY_SOLVED
    #    @test isapprox(solution["objective"], 1.36e6, rtol=1.0e-4)
    #end

    @testset "Shamir network (unknown flow directions), MILP formulation." begin
        data = WaterModels.parse_file("../test/data/epanet/shamir.inp")
        modifications = WaterModels.parse_file("../test/data/json/shamir.json")
        InfrastructureModels.update_data!(data, modifications)
        wm = build_model(data, MILPWaterModel, WaterModels.post_ne, ext=Dict(:num_breakpoints => 5))
        solution = WaterModels.optimize_model!(wm, cbc)
        @test solution["termination_status"] == OPTIMAL
        @test isapprox(solution["objective"], 1.36e6, rtol=1.0e-4)
    end

    #@testset "Shamir network (unknown flow directions), MILPR formulation." begin
    #    data = WaterModels.parse_file("../test/data/epanet/shamir.inp")
    #    modifications = WaterModels.parse_file("../test/data/json/shamir-reduced.json")
    #    InfrastructureModels.update_data!(data, modifications)
    #    wm = build_model(data, MILPRWaterModel, WaterModels.post_ne, ext=Dict(:num_breakpoints => 5))
    #    solution = WaterModels.optimize_model!(wm, cbc)
    #    @test solution["termination_status"] == OPTIMAL
    #    @test isapprox(solution["objective"], 1.36e6, rtol=1.0e-4)
    #end

    #@testset "Shamir network (reduced), NCNLP formulation." begin
    #    data = WaterModels.parse_file("../test/data/epanet/shamir.inp")
    #    modifications = WaterModels.parse_file("../test/data/json/shamir-reduced.json")
    #    InfrastructureModels.update_data!(data, modifications)
    #    wm = build_model(data, NCNLPWaterModel, WaterModels.post_ne)
    #    f = Juniper.register(fun(wm, :head_loss)..., autodiff=false)
    #    juniper = JuMP.with_optimizer(Juniper.Optimizer, nl_solver=ipopt, registered_functions=[f], log_levels=[])
    #    solution = WaterModels.optimize_model!(wm, juniper)

    #    @test solution["termination_status"] == LOCALLY_SOLVED
    #    @test isapprox(solution["objective"], 1.36e6, rtol=1.0e-4)
    #end
end
