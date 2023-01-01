@testset "src/core/node.jl" begin
    @testset "_relax_nodes! (single network with `time_series`)" begin
        data = WaterModels.parse_file("../test/data/epanet/multinetwork/reservoir-hw-lps.inp")
        WaterModels._relax_nodes!(data)
        @test isapprox(data["node"]["1"]["head_min"] * data["base_head"], 10.0)
        @test isapprox(data["node"]["1"]["head_max"] * data["base_head"], 20.0)
    end
end


@testset "make_per_unit! (node)" begin
    @testset "node, single network" begin
        # Load the data without making a per-unit transformation.
        network_path = "../test/data/epanet/multinetwork/reservoir-hw-lps.inp"
        si_data = WaterModels.parse_file(network_path; per_unit=false)

        # Transform the data to a per-unit system.
        pu_data = deepcopy(si_data)
        make_per_unit!(pu_data)

        # Ensure select pipe quantities are correctly transformed.
        head_min_in_si = pu_data["node"]["1"]["head_min"] * pu_data["base_head"]
        @test isapprox(head_min_in_si, si_data["node"]["1"]["head_min"])
        head_max_in_si = pu_data["node"]["1"]["head_max"] * pu_data["base_head"]
        @test isapprox(head_max_in_si, si_data["node"]["1"]["head_max"])
        head_nominal_in_si = pu_data["node"]["1"]["head_nominal"] * pu_data["base_head"]
        @test isapprox(head_nominal_in_si, si_data["node"]["1"]["head_nominal"])
        elevation_in_si = pu_data["node"]["1"]["elevation"] * pu_data["base_head"]
        @test isapprox(elevation_in_si, si_data["node"]["1"]["elevation"])
        head_min_ts_in_si = pu_data["time_series"]["node"]["1"]["head_min"] .* pu_data["base_head"]
        @test isapprox(head_min_ts_in_si, si_data["time_series"]["node"]["1"]["head_min"])
        head_max_ts_in_si = pu_data["time_series"]["node"]["1"]["head_max"] .* pu_data["base_head"]
        @test isapprox(head_max_ts_in_si, si_data["time_series"]["node"]["1"]["head_max"])
        head_nominal_ts_in_si = pu_data["time_series"]["node"]["1"]["head_nominal"] .* pu_data["base_head"]
        @test isapprox(head_nominal_ts_in_si, si_data["time_series"]["node"]["1"]["head_nominal"])
    end

    @testset "node, multinetwork" begin
        # Load the data without making a per-unit transformation.
        network_path = "../test/data/epanet/multinetwork/reservoir-hw-lps.inp"
        si_data = WaterModels.parse_file(network_path; per_unit=false)
        si_data_mn = make_multinetwork(si_data)

        # Transform the data to a per-unit system.
        pu_data_mn = deepcopy(si_data_mn)
        make_per_unit!(pu_data_mn)

        # Ensure select node quantities are correctly transformed.
        head_min_in_si = pu_data_mn["nw"]["1"]["node"]["1"]["head_min"] * pu_data_mn["base_head"]
        @test isapprox(head_min_in_si, si_data_mn["nw"]["1"]["node"]["1"]["head_min"])
        head_max_in_si = pu_data_mn["nw"]["1"]["node"]["1"]["head_max"] * pu_data_mn["base_head"]
        @test isapprox(head_max_in_si, si_data_mn["nw"]["1"]["node"]["1"]["head_max"])
        head_nominal_in_si = pu_data_mn["nw"]["1"]["node"]["1"]["head_nominal"] * pu_data_mn["base_head"]
        @test isapprox(head_nominal_in_si, si_data_mn["nw"]["1"]["node"]["1"]["head_nominal"])
        elevation_in_si = pu_data_mn["nw"]["1"]["node"]["1"]["elevation"] * pu_data_mn["base_head"]
        @test isapprox(elevation_in_si, si_data_mn["nw"]["1"]["node"]["1"]["elevation"])
    end
end


@testset "make_si_units! (node)" begin
    @testset "node, single network" begin
        # Load the data without making a per-unit transformation.
        network_path = "../test/data/epanet/multinetwork/reservoir-hw-lps.inp"
        si_data = WaterModels.parse_file(network_path; per_unit=false)

        # Transform the data to a per-unit system.
        data = deepcopy(si_data)
        make_per_unit!(data)

        # Transform per-unit data back to SI units.
        make_si_units!(data)

        # Ensure select node quantities are correctly transformed.
        @test isapprox(data["node"]["1"]["head_min"], si_data["node"]["1"]["head_min"])
        @test isapprox(data["node"]["1"]["head_max"], si_data["node"]["1"]["head_max"])
        @test isapprox(data["node"]["1"]["head_nominal"], si_data["node"]["1"]["head_nominal"])
        @test isapprox(data["node"]["1"]["elevation"], si_data["node"]["1"]["elevation"])
        @test isapprox(data["time_series"]["node"]["1"]["head_min"], si_data["time_series"]["node"]["1"]["head_min"])
        @test isapprox(data["time_series"]["node"]["1"]["head_max"], si_data["time_series"]["node"]["1"]["head_max"])
        @test isapprox(data["time_series"]["node"]["1"]["head_nominal"], si_data["time_series"]["node"]["1"]["head_nominal"])
    end

    @testset "node, multinetwork" begin
        # Load the data without making a per-unit transformation.
        network_path = "../test/data/epanet/multinetwork/reservoir-hw-lps.inp"
        si_data = WaterModels.parse_file(network_path; per_unit=false)
        si_data_mn = make_multinetwork(si_data)

        # Transform the data to a per-unit system.
        data_mn = deepcopy(si_data_mn)
        make_per_unit!(data_mn)

        # Transform per-unit data back to SI units.
        make_si_units!(data_mn)

        # Ensure select pipe quantities are correctly transformed.
        @test isapprox(data_mn["nw"]["1"]["node"]["1"]["head_min"], si_data_mn["nw"]["1"]["node"]["1"]["head_min"])
        @test isapprox(data_mn["nw"]["1"]["node"]["1"]["head_max"], si_data_mn["nw"]["1"]["node"]["1"]["head_max"])
        @test isapprox(data_mn["nw"]["1"]["node"]["1"]["head_nominal"], si_data_mn["nw"]["1"]["node"]["1"]["head_nominal"])
        @test isapprox(data_mn["nw"]["1"]["node"]["1"]["elevation"], si_data_mn["nw"]["1"]["node"]["1"]["elevation"])
    end
end