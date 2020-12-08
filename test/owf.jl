# Iterate over all possible WaterModels formulations.
for formulation in [NCWaterModel, CRDWaterModel, LAWaterModel, LRDWaterModel, QRDWaterModel, CQRDWaterModel]
    # Set a generic extensions dictionary for testing purposes.
    ext = Dict(:pipe_breakpoints => 3, :pump_breakpoints => 3)

    @testset "Optimal Water Flow Problems (Single Network): $(formulation)" begin
        network = WaterModels.parse_file("../test/data/epanet/snapshot/pump-hw-lps.inp")
        wm = instantiate_model(network, formulation, build_owf; ext = ext)
        result = WaterModels.optimize_model!(wm, optimizer = _make_juniper(wm, ipopt))

        @test result["termination_status"] == LOCALLY_SOLVED
        @test isapprox(result["solution"]["node"]["1"]["h"], 10.0, rtol = 1.0e-3)
        @test isapprox(result["solution"]["node"]["2"]["h"], 98.98, rtol = 2.5e-1)
        @test isapprox(result["solution"]["pump"]["1"]["status"], 1.0, atol = 1.0e-3)
        @test isapprox(result["objective"], 128.302, rtol = 5.0e-1)
    end

    @testset "Optimal Water Flow Problems (Multinetwork): $(formulation)" begin
        network = WaterModels.parse_file("../test/data/epanet/multinetwork/owf-hw-lps.inp")
        network = WaterModels.make_multinetwork(network)
        wm = instantiate_model(network, formulation, build_mn_owf; ext = ext)
        result = WaterModels.optimize_model!(wm, optimizer = _make_juniper(wm, ipopt))
        
        @test result["termination_status"] == LOCALLY_SOLVED
    end
end 