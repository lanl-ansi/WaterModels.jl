@testset "src/io/common.jl" begin
    @testset "parse_file (.inp)" begin
        data = WaterModels.parse_file("../test/data/epanet/shamir.inp")
        title = "shamir -- Bragalli, D'Ambrosio, Lee, Lodi, Toth (2008)"
        @test data["title"] == lowercase(title)
    end

    @testset "parse_file (.json)" begin
        data = WaterModels.parse_file("../test/data/json/shamir.json")
        pipe_1_max_velocity = data["pipes"]["1"]["maximumVelocity"]
        @test pipe_1_max_velocity == 2.0
    end

    @testset "parse_file (invalid extension)" begin
        path = "../test/data/json/shamir.data"
        @test_throws ErrorException WaterModels.parse_file(path)
    end

    @testset "parse_json" begin
        data = WaterModels.parse_json("../test/data/json/shamir.json")
        pipe_1_max_velocity = data["pipes"]["1"]["maximumVelocity"]
        @test pipe_1_max_velocity == 2.0
    end

    #@testset "print_solution" begin
    #    pipes = Dict{String, Any}("1" => Dict{String, Float64}("q" => 1.0))
    #    solution = Dict{String, Any}("multinetwork" => false, "pipes" => pipes)
    #    WaterModels.print_solution(solution)
    #end
end
