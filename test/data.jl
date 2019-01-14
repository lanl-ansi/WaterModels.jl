@testset "test data summary" begin
    @testset "shamir from EPANET file" begin
        data = WaterModels.parse_file("../test/data/epanet/shamir.inp")
        title = "shamir -- Bragalli, D'Ambrosio, Lee, Lodi, Toth (2008)"
        @test data["title"] == lowercase(title)
    end
end
