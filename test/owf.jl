@testset "Optimal Water Flow Problems (Single Network)" begin
    network = WaterModels.parse_file("../test/data/epanet/snapshot/pump-hw-lps.inp")

    wm = instantiate_model(network, NCWaterModel, build_owf)
    result = WaterModels.optimize_model!(wm, optimizer=_make_juniper(wm, ipopt))
    @test result["termination_status"] == LOCALLY_SOLVED
    @test isapprox(result["solution"]["node"]["1"]["h"], 10.0, rtol=1.0e-3)
    @test isapprox(result["solution"]["node"]["2"]["h"], 98.98, rtol=1.0e-3)
    @test isapprox(result["solution"]["pump"]["1"]["status"], 1.0, atol=1.0e-3)
    @test result["objective"] <= 128.302

    wm = instantiate_model(network, CRDWaterModel, build_owf, ext=Dict(:pump_breakpoints=>3))
    result = WaterModels.optimize_model!(wm, optimizer=_make_juniper(wm, ipopt))
    @test result["termination_status"] == LOCALLY_SOLVED
    @test isapprox(result["solution"]["node"]["1"]["h"], 10.0, rtol=1.0e-3)
    @test result["solution"]["node"]["2"]["h"] <= 98.99
    @test isapprox(result["solution"]["pump"]["1"]["status"], 1.0, atol=1.0e-3)
    @test result["objective"] <= 128.302

    result = run_owf(network, LAWaterModel, cbc, ext=Dict(:pump_breakpoints=>4))
    @test result["termination_status"] == OPTIMAL
    @test isapprox(result["solution"]["node"]["1"]["h"], 10.0, rtol=1.0e-3)
    @test isapprox(result["solution"]["node"]["2"]["h"], 98.98, rtol=1.0e-1)
    @test isapprox(result["solution"]["pump"]["1"]["status"], 1.0, atol=1.0e-3)
    @test result["objective"] <= 128.302

    result = run_owf(network, LRDWaterModel, cbc, ext=Dict(:pump_breakpoints=>3))
    @test result["termination_status"] == OPTIMAL
    @test isapprox(result["solution"]["node"]["1"]["h"], 10.0, rtol=1.0e-3)
    @test isapprox(result["solution"]["pump"]["1"]["status"], 1.0, atol=1.0e-3)
    @test result["objective"] <= 128.302

    wm = instantiate_model(network, QRDWaterModel, build_owf)
    result = WaterModels.optimize_model!(wm, optimizer=_make_juniper(wm, ipopt; register=false))
    @test result["termination_status"] == LOCALLY_SOLVED
    @test isapprox(result["solution"]["node"]["1"]["h"], 10.0, rtol=1.0e-3)
    @test isapprox(result["solution"]["node"]["2"]["h"], 98.98, rtol=1.0e-3)
    @test isapprox(result["solution"]["pump"]["1"]["status"], 1.0, atol=1.0e-3)

    wm = instantiate_model(network, CQRDWaterModel, build_owf)
    result = WaterModels.optimize_model!(wm, optimizer=_make_juniper(wm, ipopt; register=false))
    @test result["termination_status"] == LOCALLY_SOLVED
    @test isapprox(result["solution"]["node"]["1"]["h"], 10.0, rtol=1.0e-3)
    @test result["solution"]["node"]["2"]["h"] <= 98.99
    @test isapprox(result["solution"]["pump"]["1"]["status"], 1.0, atol=1.0e-3)
end

@testset "Optimal Water Flow Problems (Multinetwork)" begin
    network = WaterModels.parse_file("../test/data/epanet/multinetwork/owf-hw-lps.inp")
    network = WaterModels.make_multinetwork(network)

    wm = instantiate_model(network, NCWaterModel, build_mn_owf)
    result = WaterModels.optimize_model!(wm, optimizer=_make_juniper(wm, ipopt))
    # @test result["termination_status"] == LOCALLY_SOLVED # TODO: Why does this fail?

    wm = instantiate_model(network, CRDWaterModel, build_mn_owf, ext=Dict(:pump_breakpoints=>3))
    result = WaterModels.optimize_model!(wm, optimizer=_make_juniper(wm, ipopt))
    @test result["termination_status"] == LOCALLY_SOLVED

    result = run_mn_owf(network, LAWaterModel, cbc, ext=Dict(:pump_breakpoints=>4))
    @test result["termination_status"] == OPTIMAL

    result = run_mn_owf(network, LRDWaterModel, cbc, ext=Dict(:pump_breakpoints=>3))
    @test result["termination_status"] == OPTIMAL

    wm = instantiate_model(network, QRDWaterModel, build_mn_owf)
    result = WaterModels.optimize_model!(wm, optimizer=_make_juniper(wm, ipopt; register=false))
    @test result["termination_status"] == LOCALLY_SOLVED

    wm = instantiate_model(network, CQRDWaterModel, build_mn_owf)
    result = WaterModels.optimize_model!(wm, optimizer=_make_juniper(wm, ipopt; register=false))
    @test result["termination_status"] == LOCALLY_SOLVED
end
