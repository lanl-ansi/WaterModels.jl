@testset "Optimization-based Bound Tightening (Reduced Network)" begin
    network = WaterModels.parse_file("../test/data/epanet/multinetwork/pump-hw-lps.inp")
    run_obbt_owf!(network, ipopt)
    @test haskey(network["pump"]["1"], "q_min")
    @test haskey(network["pump"]["1"], "q_max")
end

@testset "Optimization-based Bound Tightening (Multinetwork)" begin
    network = WaterModels.parse_file("../test/data/epanet/multinetwork/pump-hw-lps.inp")
    network_mn = WaterModels.make_multinetwork(network)
    run_obbt_owf!(network_mn, ipopt, use_reduced_network=false)
    @test haskey(network_mn["nw"]["1"]["pump"]["1"], "q_min")
    @test haskey(network_mn["nw"]["1"]["pump"]["1"], "q_max")
end

@testset "Unbinarize (Multinetwork)" begin
    network = WaterModels.parse_file("../test/data/epanet/multinetwork/pump-hw-lps.inp")
    network_mn = WaterModels.make_multinetwork(network)
    wm = instantiate_model(network_mn, NCWaterModel, build_mn_wf)
    unbinarize_mn!(wm, 3) # Transform binary variables to continuous.
    @test !JuMP.is_binary(_IM.var(wm, 1, :z_pump, 1))
    @test !JuMP.is_binary(_IM.var(wm, 2, :z_pump, 1))
    @test !JuMP.is_binary(_IM.var(wm, 3, :z_pump, 1))
end
