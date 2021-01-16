@testset "src/io/epanet.jl" begin
    @testset "_initialize_epanet_dictionary" begin
        data = WaterModels._initialize_epanet_dictionary()
        @test haskey(data, "section") == true
        @test haskey(data["section"], "[TITLE]") == true
        @test haskey(data, "time_series") == true
        @test haskey(data, "top_comments") == true
    end

    @testset "_read_file_as_string (.inp)" begin
        file_path = "../test/data/epanet/snapshot/pipe-hw-lps.inp"
        contents = WaterModels._read_file_as_string(file_path)
        @test occursin("[TITLE]", contents) == true
        @test occursin("[PIPES]", contents) == true
    end

    @testset "_read_file_as_string (does not exist)" begin
        file_path = "/this/file/probably/does/not/exist.inp"
        @test_throws ErrorException WaterModels._read_file_as_string(file_path)
    end

    @testset "_read_epanet_sections" begin
        file_path = "../test/data/epanet/snapshot/pipe-hw-lps.inp"
        raw_data = WaterModels._read_epanet_sections(file_path)
        @test haskey(raw_data["section"], "[TITLE]") == true
        @test haskey(raw_data["section"], "[PIPES]") == true
    end

    @testset "_read_epanet_sections (top comments)" begin
        file_path = "../test/data/epanet/snapshot/top-comment.inp"
        raw_data = WaterModels._read_epanet_sections(file_path)
        @test length(raw_data["top_comments"]) > 0
    end

    @testset "_read_epanet_sections (section title error)" begin
        file_path = "../test/data/epanet/snapshot/error-section.inp"
        @test_throws ErrorException WaterModels._read_epanet_sections(file_path)
    end

    @testset "_read_epanet_sections (invalid syntax)" begin
        file_path = "../test/data/epanet/snapshot/no-sections.inp"
        @test_throws ErrorException WaterModels._read_epanet_sections(file_path)
    end
end
