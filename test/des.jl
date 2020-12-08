@testset "Network Design Problems (Single Network)" begin
    network = parse_file("../examples/data/epanet/shamir.inp")
    modifications = parse_json("../test/data/json/shamir-reduced.json")
    _IM.update_data!(network, modifications)

    @testset "Shamir network (reduced), NC formulation." begin
        wm = instantiate_model(network, NCWaterModel, build_des)
        result = WaterModels.optimize_model!(wm, optimizer=_make_juniper(wm, ipopt))
        @test result["termination_status"] == LOCALLY_SOLVED
        @test isapprox(result["objective"], 1.36e6, rtol=1.0e-4)
    end

    @testset "Shamir network (reduced), CRD formulation." begin
        wm = instantiate_model(network, CRDWaterModel, build_des)
        result = WaterModels.optimize_model!(wm, optimizer=_make_juniper(wm, ipopt))
        @test result["termination_status"] == LOCALLY_SOLVED
        @test isapprox(result["objective"], 1.36e6, rtol=1.0e-4)
    end

    @testset "Shamir network (reduced), LA formulation." begin
        result = run_des(network, LAWaterModel, cbc, ext=Dict(:pipe_breakpoints=>3))
        @test result["termination_status"] == OPTIMAL
        @test isapprox(result["objective"], 1.36e6, rtol=1.0e-4)
    end

    @testset "Shamir network (reduced), LRD formulation." begin
        result = run_des(network, LRDWaterModel, cbc, ext=Dict(:pipe_breakpoints=>3))
        @test result["termination_status"] == OPTIMAL
        @test isapprox(result["objective"], 1.36e6, rtol=1.0e-4)
    end
end
