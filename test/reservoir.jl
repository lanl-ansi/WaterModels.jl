@testset "src/core/reservoir.jl" begin
    @testset "_fix_reservoirs! (multinetwork)" begin
        data = WaterModels.parse_file("../test/data/epanet/multinetwork/pipe-hw-lps.inp")
        mn_data = WaterModels.make_multinetwork(data)

        mn_data["nw"]["1"]["reservoir"]["1"]["dispatchable"] = true
        mn_data["nw"]["3"]["reservoir"]["1"]["dispatchable"] = true

        WaterModels._fix_reservoirs!(mn_data)
        @test mn_data["nw"]["1"]["reservoir"]["1"]["dispatchable"] == false
        @test mn_data["nw"]["3"]["reservoir"]["1"]["dispatchable"] == false
    end

    @testset "_relax_reservoirs! (single network)" begin
        data = WaterModels.parse_file("../test/data/epanet/snapshot/pipe-hw-lps.inp")
        WaterModels._relax_reservoirs!(data)
        @test data["reservoir"]["1"]["dispatchable"] == true
    end

    @testset "_relax_reservoirs! (multinetwork)" begin
        data = WaterModels.parse_file("../test/data/epanet/multinetwork/pipe-hw-lps.inp")
        mn_data = WaterModels.make_multinetwork(data)
        WaterModels._relax_reservoirs!(mn_data)
        @test mn_data["nw"]["1"]["reservoir"]["1"]["dispatchable"] == true 
        @test mn_data["nw"]["3"]["reservoir"]["1"]["dispatchable"] == true
    end
end


@testset "make_per_unit! (reservoir)" begin
    @testset "reservoir, single network" begin
        # Load the data without making a per-unit transformation.
        network_path = "../test/data/epanet/multinetwork/reservoir-hw-lps.inp"
        si_data = WaterModels.parse_file(network_path; per_unit=false)

        # Transform the data to a per-unit system.
        pu_data = deepcopy(si_data)
        make_per_unit!(pu_data)

        # Ensure select reservoir quantities are correctly transformed.
        head_min_in_si = pu_data["reservoir"]["1"]["head_min"] * pu_data["base_head"]
        @test isapprox(head_min_in_si, si_data["reservoir"]["1"]["head_min"])
        head_max_in_si = pu_data["reservoir"]["1"]["head_max"] * pu_data["base_head"]
        @test isapprox(head_max_in_si, si_data["reservoir"]["1"]["head_max"])
        head_nominal_in_si = pu_data["reservoir"]["1"]["head_nominal"] * pu_data["base_head"]
        @test isapprox(head_nominal_in_si, si_data["reservoir"]["1"]["head_nominal"])
    end

    @testset "reservoir, multinetwork" begin
        # Load the data without making a per-unit transformation.
        network_path = "../test/data/epanet/multinetwork/reservoir-hw-lps.inp"
        si_data = WaterModels.parse_file(network_path; per_unit=false)
        si_data_mn = make_multinetwork(si_data)

        # Transform the data to a per-unit system.
        pu_data_mn = deepcopy(si_data_mn)
        make_per_unit!(pu_data_mn)

        # Ensure select reservoir quantities are correctly transformed.
        head_min_in_si = pu_data_mn["nw"]["1"]["reservoir"]["1"]["head_min"] * pu_data_mn["base_head"]
        @test isapprox(head_min_in_si, si_data_mn["nw"]["1"]["reservoir"]["1"]["head_min"])
        head_max_in_si = pu_data_mn["nw"]["1"]["reservoir"]["1"]["head_max"] * pu_data_mn["base_head"]
        @test isapprox(head_max_in_si, si_data_mn["nw"]["1"]["reservoir"]["1"]["head_max"])
        head_nominal_in_si = pu_data_mn["nw"]["1"]["reservoir"]["1"]["head_nominal"] * pu_data_mn["base_head"]
        @test isapprox(head_nominal_in_si, si_data_mn["nw"]["1"]["reservoir"]["1"]["head_nominal"])
    end
end


@testset "make_si_units! (reservoir)" begin
    @testset "reservoir, single network" begin
        # Load the data without making a per-unit transformation.
        network_path = "../test/data/epanet/multinetwork/reservoir-hw-lps.inp"
        si_data = WaterModels.parse_file(network_path; per_unit=false)

        # Transform the data to a per-unit system.
        data = deepcopy(si_data)
        make_per_unit!(data)

        # Transform per-unit data back to SI units.
        make_si_units!(data)

        # Ensure select reservoir quantities are correctly transformed.
        @test isapprox(data["reservoir"]["1"]["head_min"], si_data["reservoir"]["1"]["head_min"])
        @test isapprox(data["reservoir"]["1"]["head_max"], si_data["reservoir"]["1"]["head_max"])
        @test isapprox(data["reservoir"]["1"]["head_nominal"], si_data["reservoir"]["1"]["head_nominal"])

    end

    @testset "reservoir, multinetwork" begin
        # Load the data without making a per-unit transformation.
        network_path = "../test/data/epanet/multinetwork/reservoir-hw-lps.inp"
        si_data = WaterModels.parse_file(network_path; per_unit=false)
        si_data_mn = make_multinetwork(si_data)

        # Transform the data to a per-unit system.
        data_mn = deepcopy(si_data_mn)
        make_per_unit!(data_mn)

        # Transform per-unit data back to SI units.
        make_si_units!(data_mn)

        # Ensure select reservoir quantities are correctly transformed.
        @test isapprox(data_mn["nw"]["1"]["reservoir"]["1"]["head_min"], si_data_mn["nw"]["1"]["reservoir"]["1"]["head_min"])
        @test isapprox(data_mn["nw"]["1"]["reservoir"]["1"]["head_max"], si_data_mn["nw"]["1"]["reservoir"]["1"]["head_max"])
        @test isapprox(data_mn["nw"]["1"]["reservoir"]["1"]["head_nominal"], si_data_mn["nw"]["1"]["reservoir"]["1"]["head_nominal"])
    end
end
