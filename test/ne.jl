@testset "Network Expansion Problems" begin
    # Get an initial solution for the Hazen-Williams network.
    hw_network_path = "../test/data/epanet/shamir.inp"
    hw_network = WaterModels.parse_file(hw_network_path)
    hw_wf_result = WaterModels.run_wf(hw_network, CNLPWaterModel, ipopt, alpha=1.852)
    InfrastructureModels.update_data!(hw_network, hw_wf_result["solution"])

    @testset "Shamir network (unknown flow directions), CNLP formulation." begin
        network = deepcopy(hw_network)
        modifications = WaterModels.parse_file("../test/data/json/shamir.json")
        InfrastructureModels.update_data!(network, modifications)
        @test_throws ErrorException run_ne(network, CNLPWaterModel, ipopt, alpha=1.852, relaxed=true)
    end

    @testset "Shamir network (unknown flow directions), MICP formulation." begin
        network = deepcopy(hw_network)
        modifications = WaterModels.parse_file("../test/data/json/shamir.json")
        InfrastructureModels.update_data!(network, modifications)
        result = run_ne(network, MICPWaterModel, ipopt, alpha=1.852, relaxed=true)
        @test result["termination_status"] == MOI.ALMOST_LOCALLY_SOLVED ||
              result["termination_status"] == MOI.LOCALLY_SOLVED
    end

    #@testset "Shamir network (unknown flow directions), NCNLP formulation." begin
    #    network = WaterModels.parse_file("../test/data/epanet/shamir.inp")
    #    modifications = WaterModels.parse_file("../test/data/json/shamir.json")
    #    InfrastructureModels.update_data!(network, modifications)
    #    solution = run_ne(network, NCNLPWaterModel, ipopt, alpha=1.852, relaxed=true)
    #    @test solution["termination_status"] == MOI.ALMOST_LOCALLY_SOLVED ||
    #          solution["termination_status"] == MOI.LOCALLY_SOLVED
    #end
end
