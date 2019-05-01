@testset "Network Expansion Problems" begin
    #@testset "Shamir network (unknown flow directions), MICP formulation." begin
    #    network = WaterModels.parse_file("../test/data/epanet/shamir.inp")
    #    modifications = WaterModels.parse_file("../test/data/json/shamir.json")
    #    InfrastructureModels.update_data!(network, modifications)
    #    solution = run_ne(network, MICPWaterModel, ipopt, alpha=1.852, relaxed=true)
    #    @test solution["termination_status"] == MOI.ALMOST_LOCALLY_SOLVED ||
    #          solution["termination_status"] == MOI.LOCALLY_SOLVED
    #end

    #@testset "Shamir network (unknown flow directions), NCNLP formulation." begin
    #    network = WaterModels.parse_file("../test/data/epanet/shamir.inp")
    #    modifications = WaterModels.parse_file("../test/data/json/shamir.json")
    #    InfrastructureModels.update_data!(network, modifications)
    #    solution = run_ne(network, NCNLPWaterModel, ipopt, alpha=1.852, relaxed=true)
    #    @test solution["termination_status"] == MOI.ALMOST_LOCALLY_SOLVED ||
    #          solution["termination_status"] == MOI.LOCALLY_SOLVED
    #end
end
