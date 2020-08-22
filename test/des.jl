@testset "Network Design Problems (Single Network)" begin
    network = parse_file("../examples/data/epanet/shamir.inp")
    modifications = parse_file("../test/data/json/shamir-reduced.json")
    _IM.update_data!(network, modifications)

    @testset "Hazen-Williams NLP Formulation." begin
        wm = instantiate_model(network, NLPWaterModel, build_des)
        result = WaterModels.optimize_model!(wm, optimizer=_make_juniper(wm, ipopt))
        @test result["termination_status"] == LOCALLY_SOLVED
        @test isapprox(result["objective"], 1.36e6, rtol=1.0e-4)
    end

    @testset "Hazen-Williams NLP Formulation." begin
        wm = instantiate_model(network, MICPRWaterModel, build_des)
        result = WaterModels.optimize_model!(wm, optimizer=_make_juniper(wm, ipopt))
        @test result["termination_status"] == LOCALLY_SOLVED
        @test isapprox(result["objective"], 1.36e6, rtol=1.0e-4)
    end

    @testset "Shamir network (reduced), MILP formulation." begin
        result = run_des(network, MILPWaterModel, cbc, ext=Dict(:pipe_breakpoints=>3))
        @test result["termination_status"] == OPTIMAL
        @test isapprox(result["objective"], 1.36e6, rtol=1.0e-4)
    end

    @testset "Shamir network (reduced), MILP-R formulation." begin
        result = run_des(network, MILPRWaterModel, cbc, ext=Dict(:pipe_breakpoints=>3))
        @test result["termination_status"] == OPTIMAL
        @test isapprox(result["objective"], 1.36e6, rtol=1.0e-4)
    end
end
