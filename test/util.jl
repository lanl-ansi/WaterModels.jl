@testset "Optimization-based Bound Tightening (Reduced Network)" begin
    network = WaterModels.parse_file("../test/data/epanet/multinetwork/pump-hw-lps.inp")
    run_obbt_owf!(network, ipopt)
    @test haskey(network["pump"]["1"], "q_max")
end
