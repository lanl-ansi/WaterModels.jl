@testset "Hazen-Williams NLP Problems" begin
    @testset "Hanoi network." begin
        network_path = "../test/data/epanet/hanoi.inp"
        modification_path = "../test/data/json/ne-hanoi.json"
        solution = run_ne_hw(network_path, modification_path, NLPWaterModel, knitro)
        println(solution["objective"])
        @test solution["status"] == :LocalOptimal
    end
end
