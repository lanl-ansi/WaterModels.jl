@testset "test data summary" begin
    @testset "shamir from EPANET file" begin
        data = WaterModels.parse_file("../test/data/epanet/shamir.inp")
        @test data["title"] == lowercase("shamir -- Bragalli, D'Ambrosio, Lee, Lodi, Toth (2008)")
    end
end
