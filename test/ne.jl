# Iterate over all possible WaterModels formulations.
for formulation in [NCWaterModel, NCDWaterModel, CRDWaterModel, LAWaterModel, LRDWaterModel, PWLRDWaterModel]
    @testset "Network Expansion Problems (Single Network): $(formulation)" begin
        network = WaterModels.parse_file("../test/data/epanet/snapshot/pipe-hw-lps.inp"; skip_correct = true)
        network["ne_short_pipe"]["1"] = pop!(network["pipe"], "1")

        delete!.(Ref(network["ne_short_pipe"]["1"]), ["diameter", "length", "has_valve", "roughness"])
        network["ne_short_pipe"]["1"]["status"] = WaterModels.STATUS_UNKNOWN
        network["ne_short_pipe"]["1"]["construction_cost"] = 10.0
        WaterModels.correct_network_data!(network)

        wm = instantiate_model(network, formulation, build_ne)
        result = WaterModels.optimize_model!(wm, optimizer = _choose_solver(wm, nlp_solver, milp_solver))

        @test _is_valid_status(result["termination_status"])
        @test isapprox(result["solution"]["ne_short_pipe"]["1"]["status"], 1.0, atol = 1.0e-3)
        @test isapprox(result["objective"], 10.0)
    end

    @testset "Network Expansions Problems (Multinetwork): $(formulation)" begin
        network = WaterModels.parse_file("../test/data/epanet/multinetwork/pipe-hw-lps.inp"; skip_correct = true)
        network["ne_short_pipe"]["1"] = pop!(network["pipe"], "1")

        delete!.(Ref(network["ne_short_pipe"]["1"]), ["diameter", "length", "has_valve", "roughness"])
        network["ne_short_pipe"]["1"]["status"] = WaterModels.STATUS_UNKNOWN
        network["ne_short_pipe"]["1"]["construction_cost"] = 10.0
        WaterModels.correct_network_data!(network)
        network_mn = WaterModels.make_multinetwork(network)

        wm = instantiate_model(network_mn, formulation, build_mn_ne)
        result = WaterModels.optimize_model!(wm, optimizer = _choose_solver(wm, nlp_solver, milp_solver))

        @test _is_valid_status(result["termination_status"])
        @test isapprox(result["solution"]["nw"]["1"]["ne_short_pipe"]["1"]["status"], 1.0, atol = 1.0e-3)
        @test isapprox(result["solution"]["nw"]["2"]["ne_short_pipe"]["1"]["status"], 1.0, atol = 1.0e-3)
        @test isapprox(result["objective"], (length(result["solution"]["nw"]) - 1.0) * 10.0)
    end
end

@testset "solve_ne" begin
    network = WaterModels.parse_file("../test/data/epanet/snapshot/pipe-hw-lps.inp"; skip_correct = true)
    network["ne_short_pipe"]["1"] = pop!(network["pipe"], "1")

    delete!.(Ref(network["ne_short_pipe"]["1"]), ["diameter", "length", "has_valve", "roughness"])
    network["ne_short_pipe"]["1"]["status"] = WaterModels.STATUS_UNKNOWN
    network["ne_short_pipe"]["1"]["construction_cost"] = 10.0
    WaterModels.correct_network_data!(network)

    result = WaterModels.solve_ne(network, LRDWaterModel, milp_solver)
    result = WaterModels.run_ne(network, LRDWaterModel, milp_solver)
    @test _is_valid_status(result["termination_status"])
end


@testset "solve_mn_ne" begin
    network = WaterModels.parse_file("../test/data/epanet/multinetwork/pipe-hw-lps.inp"; skip_correct = true)
    network["ne_short_pipe"]["1"] = pop!(network["pipe"], "1")

    delete!.(Ref(network["ne_short_pipe"]["1"]), ["diameter", "length", "has_valve", "roughness"])
    network["ne_short_pipe"]["1"]["status"] = WaterModels.STATUS_UNKNOWN
    network["ne_short_pipe"]["1"]["construction_cost"] = 10.0
    WaterModels.correct_network_data!(network)

    network_mn = WaterModels.make_multinetwork(network)
    result = WaterModels.solve_mn_ne(network_mn, LRDWaterModel, milp_solver)
    result = WaterModels.run_mn_ne(network_mn, LRDWaterModel, milp_solver)
    @test _is_valid_status(result["termination_status"])
end