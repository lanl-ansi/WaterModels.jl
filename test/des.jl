@testset "Network Design Problems" begin
    @testset "Shamir network (reduced), MICP-R formulation." begin
        data = parse_file("../test/data/epanet/shamir.inp")
        modifications = parse_file("../test/data/json/shamir-reduced.json")
        InfrastructureModels.update_data!(data, modifications)
        wm = instantiate_model(data, MICPRWaterModel, build_des)

        f = Juniper.register(head_loss_args(wm)..., autodiff=false)
        juniper = JuMP.optimizer_with_attributes(Juniper.Optimizer,
            "nl_solver"=>ipopt, "registered_functions"=>[f], "log_levels"=>[])
        solution = _IM.optimize_model!(wm, optimizer=juniper)

        @test solution["termination_status"] == LOCALLY_SOLVED
        @test isapprox(solution["objective"], 1.36e6, rtol=1.0e-4)
    end

    @testset "Shamir network (reduced), MILP formulation." begin
        data = parse_file("../test/data/epanet/shamir.inp")
        modifications = parse_file("../test/data/json/shamir-reduced.json")
        InfrastructureModels.update_data!(data, modifications)

        wm = instantiate_model(data, MILPWaterModel, build_des)
        solution = _IM.optimize_model!(wm, optimizer=cbc)

        @test solution["termination_status"] == OPTIMAL
        @test isapprox(solution["objective"], 1.36e6, rtol=1.0e-4)
    end

    @testset "Shamir network (reduced), MILP-R formulation." begin
        data = parse_file("../test/data/epanet/shamir.inp")
        modifications = parse_file("../test/data/json/shamir-reduced.json")
        InfrastructureModels.update_data!(data, modifications)

        wm = instantiate_model(data, MILPRWaterModel, build_des)
        solution = _IM.optimize_model!(wm, optimizer=cbc)

        @test solution["termination_status"] == OPTIMAL
        @test isapprox(solution["objective"], 1.36e6, rtol=1.0e-4)
    end

    @testset "Shamir network (reduced), NLP formulation." begin
        data = parse_file("../test/data/epanet/shamir.inp")
        modifications = parse_file("../test/data/json/shamir-reduced.json")
        InfrastructureModels.update_data!(data, modifications)
        wm = instantiate_model(data, NLPWaterModel, build_des)

        f = Juniper.register(head_loss_args(wm)..., autodiff=false)
        juniper = JuMP.optimizer_with_attributes(Juniper.Optimizer,
            "nl_solver"=>ipopt, "registered_functions"=>[f], "log_levels"=>[])
        solution = _IM.optimize_model!(wm, optimizer=juniper)

        @test solution["termination_status"] == LOCALLY_SOLVED
        @test isapprox(solution["objective"], 1.36e6, rtol=1.0e-4)
    end
end
