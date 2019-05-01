@testset "Hazen-Williams MICP Network Expansion Problems" begin
    @testset "Shamir network (unknown flow directions)." begin
        network = WaterModels.parse_file("../test/data/epanet/shamir.inp")
        modifications = WaterModels.parse_file("../test/data/json/shamir.json")
        InfrastructureModels.update_data!(network, modifications)
        solution = run_ne(network, MICPWaterModel, ipopt, alpha=1.852, relaxed=true)
        @test solution["termination_status"] == MOI.ALMOST_LOCALLY_SOLVED ||
              solution["termination_status"] == MOI.LOCALLY_SOLVED
    end
end
