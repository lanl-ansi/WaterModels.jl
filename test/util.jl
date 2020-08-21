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
