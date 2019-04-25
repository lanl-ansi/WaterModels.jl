@testset "Hazen-Williams MILP-R Network Expansion Problems" begin
    #@testset "Shamir network (unknown flow directions)." begin
    #    network = WaterModels.parse_file("../test/data/epanet/shamir.inp")
    #    modifications = WaterModels.parse_file("../test/data/json/shamir.json")
    #    InfrastructureModels.update_data!(network, modifications)
    #    solution = run_ne(network, MILPRWaterModel, cbc, alpha=0.852)
    #    @test solution["termination_status"] == MOI.OPTIMAL
    #end
end
