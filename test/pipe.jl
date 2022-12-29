@testset "make_per_unit! (pipe)" begin
    @testset "pipe, single network" begin
        # Load the data without making a per-unit transformation.
        network_path = "../test/data/epanet/snapshot/pipe-hw-lps.inp"
        si_data = WaterModels.parse_file(network_path; per_unit=false)
        si_data["pipe"]["1"]["minor_loss"] = 1.0

        # Transform the data to a per-unit system.
        pu_data = deepcopy(si_data)
        make_per_unit!(pu_data)

        # Ensure select pipe quantities are correctly transformed.
        length_in_si = pu_data["pipe"]["1"]["length"] * pu_data["base_length"]
        @test isapprox(length_in_si, si_data["pipe"]["1"]["length"])
        diameter_in_si = pu_data["pipe"]["1"]["diameter"] * pu_data["base_length"]
        @test isapprox(diameter_in_si, si_data["pipe"]["1"]["diameter"])
        flow_min_in_si = pu_data["pipe"]["1"]["flow_min"] * pu_data["base_flow"]
        @test isapprox(flow_min_in_si, si_data["pipe"]["1"]["flow_min"])
        flow_max_in_si = pu_data["pipe"]["1"]["flow_max"] * pu_data["base_flow"]
        @test isapprox(flow_max_in_si, si_data["pipe"]["1"]["flow_max"])
        minor_loss_in_si = pu_data["pipe"]["1"]["minor_loss"] * pu_data["base_flow"]
        @test isapprox(minor_loss_in_si, si_data["pipe"]["1"]["minor_loss"])
    end

    @testset "pipe, multinetwork" begin
        # Load the data without making a per-unit transformation.
        network_path = "../test/data/epanet/multinetwork/pipe-hw-lps.inp"
        si_data = WaterModels.parse_file(network_path; per_unit=false)
        si_data_mn = make_multinetwork(si_data)
        si_data_mn["nw"]["2"]["pipe"]["1"]["minor_loss"] = 1.0

        # Transform the data to a per-unit system.
        pu_data_mn = deepcopy(si_data_mn)
        make_per_unit!(pu_data_mn)

        # Ensure select pipe quantities are correctly transformed.
        length_in_si = pu_data_mn["nw"]["2"]["pipe"]["1"]["length"] * pu_data_mn["base_length"]
        @test isapprox(length_in_si, si_data_mn["nw"]["2"]["pipe"]["1"]["length"])
        diameter_in_si = pu_data_mn["nw"]["2"]["pipe"]["1"]["diameter"] * pu_data_mn["base_length"]
        @test isapprox(diameter_in_si, si_data_mn["nw"]["2"]["pipe"]["1"]["diameter"])
        flow_min_in_si = pu_data_mn["nw"]["2"]["pipe"]["1"]["flow_min"] * pu_data_mn["base_flow"]
        @test isapprox(flow_min_in_si, si_data_mn["nw"]["2"]["pipe"]["1"]["flow_min"])
        flow_max_in_si = pu_data_mn["nw"]["2"]["pipe"]["1"]["flow_max"] * pu_data_mn["base_flow"]
        @test isapprox(flow_max_in_si, si_data_mn["nw"]["2"]["pipe"]["1"]["flow_max"])
        minor_loss_in_si = pu_data_mn["nw"]["2"]["pipe"]["1"]["minor_loss"] * pu_data_mn["base_flow"]
        @test isapprox(minor_loss_in_si, si_data_mn["nw"]["2"]["pipe"]["1"]["minor_loss"])
    end

    @testset "des_pipe, single network" begin
        # Load the data without making a per-unit transformation.
        network_path = "../examples/data/json/shamir.json"
        si_data = WaterModels.parse_file(network_path; per_unit=false)
        si_data["des_pipe"]["1"]["minor_loss"] = 1.0

        # Transform the data to a per-unit system.
        pu_data = deepcopy(si_data)
        make_per_unit!(pu_data)

        # Ensure select pipe quantities are correctly transformed.
        length_in_si = pu_data["des_pipe"]["1"]["length"] * pu_data["base_length"]
        @test isapprox(length_in_si, si_data["des_pipe"]["1"]["length"])
        diameter_in_si = pu_data["des_pipe"]["1"]["diameter"] * pu_data["base_length"]
        @test isapprox(diameter_in_si, si_data["des_pipe"]["1"]["diameter"])
        flow_min_in_si = pu_data["des_pipe"]["1"]["flow_min"] * pu_data["base_flow"]
        @test isapprox(flow_min_in_si, si_data["des_pipe"]["1"]["flow_min"])
        flow_max_in_si = pu_data["des_pipe"]["1"]["flow_max"] * pu_data["base_flow"]
        @test isapprox(flow_max_in_si, si_data["des_pipe"]["1"]["flow_max"])
        minor_loss_in_si = pu_data["des_pipe"]["1"]["minor_loss"] * pu_data["base_flow"]
        @test isapprox(minor_loss_in_si, si_data["des_pipe"]["1"]["minor_loss"])
    end
end

@testset "make_si_units! (pipe)" begin
    @testset "pipe, single network" begin
        # Load the data without making a per-unit transformation.
        network_path = "../test/data/epanet/snapshot/pipe-hw-lps.inp"
        si_data = WaterModels.parse_file(network_path; per_unit=false)
        si_data["pipe"]["1"]["minor_loss"] = 1.0

        # Transform the data to a per-unit system.
        data = deepcopy(si_data)
        make_per_unit!(data)

        # Transform per-unit data back to SI units.
        make_si_units!(data)

        # Ensure select pipe quantities are correctly transformed.
        @test isapprox(data["pipe"]["1"]["length"], si_data["pipe"]["1"]["length"])
        @test isapprox(data["pipe"]["1"]["diameter"], si_data["pipe"]["1"]["diameter"])
        @test isapprox(data["pipe"]["1"]["flow_min"], si_data["pipe"]["1"]["flow_min"])
        @test isapprox(data["pipe"]["1"]["flow_max"], si_data["pipe"]["1"]["flow_max"])
        @test isapprox(data["pipe"]["1"]["minor_loss"], si_data["pipe"]["1"]["minor_loss"])
    end

    @testset "pipe, multinetwork" begin
        # Load the data without making a per-unit transformation.
        network_path = "../test/data/epanet/multinetwork/pipe-hw-lps.inp"
        si_data = WaterModels.parse_file(network_path; per_unit=false)
        si_data_mn = make_multinetwork(si_data)
        si_data_mn["nw"]["2"]["pipe"]["1"]["minor_loss"] = 1.0

        # Transform the data to a per-unit system.
        data_mn = deepcopy(si_data_mn)
        make_per_unit!(data_mn)

        # Transform per-unit data back to SI units.
        make_si_units!(data_mn)

        # Ensure select pipe quantities are correctly transformed.
        @test isapprox(data_mn["nw"]["2"]["pipe"]["1"]["length"], si_data_mn["nw"]["2"]["pipe"]["1"]["length"])
        @test isapprox(data_mn["nw"]["2"]["pipe"]["1"]["diameter"], si_data_mn["nw"]["2"]["pipe"]["1"]["diameter"])
        @test isapprox(data_mn["nw"]["2"]["pipe"]["1"]["flow_min"], si_data_mn["nw"]["2"]["pipe"]["1"]["flow_min"])
        @test isapprox(data_mn["nw"]["2"]["pipe"]["1"]["flow_max"], si_data_mn["nw"]["2"]["pipe"]["1"]["flow_max"])
        @test isapprox(data_mn["nw"]["2"]["pipe"]["1"]["minor_loss"], si_data_mn["nw"]["2"]["pipe"]["1"]["minor_loss"])
    end

    @testset "des_pipe, single network" begin
        # Load the data without making a per-unit transformation.
        network_path = "../examples/data/json/shamir.json"
        si_data = WaterModels.parse_file(network_path; per_unit=false)
        si_data["des_pipe"]["1"]["minor_loss"] = 1.0

        # Transform the data to a per-unit system.
        data = deepcopy(si_data)
        make_per_unit!(data)

        # Transform per-unit data back to SI units.
        make_si_units!(data)

        # Ensure select pipe quantities are correctly transformed.
        @test isapprox(data["des_pipe"]["1"]["length"], si_data["des_pipe"]["1"]["length"])
        @test isapprox(data["des_pipe"]["1"]["diameter"], si_data["des_pipe"]["1"]["diameter"])
        @test isapprox(data["des_pipe"]["1"]["flow_min"], si_data["des_pipe"]["1"]["flow_min"])
        @test isapprox(data["des_pipe"]["1"]["flow_max"], si_data["des_pipe"]["1"]["flow_max"])
        @test isapprox(data["des_pipe"]["1"]["minor_loss"], si_data["des_pipe"]["1"]["minor_loss"])
    end
end