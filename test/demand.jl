@testset "src/core/demand.jl" begin
    @testset "_fix_demands! (multinetwork)" begin
        data = WaterModels.parse_file("../test/data/epanet/multinetwork/pipe-hw-lps.inp")
        mn_data = WaterModels.make_multinetwork(data)
        WaterModels._fix_demands!(mn_data)
        @test mn_data["nw"]["1"]["demand"]["2"]["dispatchable"] == false
        @test mn_data["nw"]["3"]["demand"]["2"]["dispatchable"] == false
    end

    @testset "_relax_demand! (has keys)" begin
        data = WaterModels.parse_file("../test/data/epanet/snapshot/pipe-hw-lps.inp")
        data["demand"]["2"]["flow_min"], data["demand"]["2"]["flow_max"] = 0.0, 1.0
        WaterModels._relax_demand!(data["demand"]["2"])
        @test data["demand"]["2"]["dispatchable"] == true
    end

    @testset "_relax_demand! (does not have keys)" begin
        data = WaterModels.parse_file("../test/data/epanet/snapshot/pipe-hw-lps.inp")
        WaterModels._relax_demand!(data["demand"]["2"])
        @test data["demand"]["2"]["dispatchable"] == false
    end

    @testset "_relax_demands! (single network)" begin
        data = WaterModels.parse_file("../test/data/epanet/snapshot/pipe-hw-lps.inp")
        WaterModels._relax_demands!(data)
        @test data["demand"]["2"]["dispatchable"] == false
    end

    @testset "_relax_demands! (single network with `time_series`)" begin
        data = WaterModels.parse_file("../test/data/epanet/multinetwork/pipe-hw-lps.inp")
        WaterModels._relax_demands!(data)
        @test data["demand"]["2"]["dispatchable"] == true
    end

    @testset "_relax_demands! (multinetwork)" begin
        data = WaterModels.parse_file("../test/data/epanet/multinetwork/pipe-hw-lps.inp")
        mn_data = WaterModels.make_multinetwork(data)
        WaterModels._relax_demands!(mn_data)
        @test mn_data["nw"]["1"]["demand"]["2"]["dispatchable"] == false
        @test mn_data["nw"]["3"]["demand"]["2"]["dispatchable"] == false
    end
end


@testset "make_per_unit! (demand)" begin
    @testset "demand, single network" begin
        # Load the data without making a per-unit transformation.
        network_path = "../test/data/epanet/multinetwork/reservoir-hw-lps.inp"
        si_data = WaterModels.parse_file(network_path; per_unit=false)

        # Transform the data to a per-unit system.
        pu_data = deepcopy(si_data)
        make_per_unit!(pu_data)

        # Ensure select demand quantities are correctly transformed.
        flow_min_in_si = pu_data["demand"]["2"]["flow_min"] * pu_data["base_flow"]
        @test isapprox(flow_min_in_si, si_data["demand"]["2"]["flow_min"])
        flow_max_in_si = pu_data["demand"]["2"]["flow_max"] * pu_data["base_flow"]
        @test isapprox(flow_max_in_si, si_data["demand"]["2"]["flow_max"])
        flow_nominal_in_si = pu_data["demand"]["2"]["flow_nominal"] * pu_data["base_flow"]
        @test isapprox(flow_nominal_in_si, si_data["demand"]["2"]["flow_nominal"])
        flow_min_ts_in_si = pu_data["time_series"]["demand"]["2"]["flow_min"] * pu_data["base_flow"]
        @test isapprox(flow_min_ts_in_si, si_data["time_series"]["demand"]["2"]["flow_min"])
        flow_max_ts_in_si = pu_data["time_series"]["demand"]["2"]["flow_max"] * pu_data["base_flow"]
        @test isapprox(flow_max_ts_in_si, si_data["time_series"]["demand"]["2"]["flow_max"])
        flow_nominal_ts_in_si = pu_data["time_series"]["demand"]["2"]["flow_nominal"] * pu_data["base_flow"]
        @test isapprox(flow_nominal_ts_in_si, si_data["time_series"]["demand"]["2"]["flow_nominal"])
    end

    @testset "demand, multinetwork" begin
        # Load the data without making a per-unit transformation.
        network_path = "../test/data/epanet/multinetwork/reservoir-hw-lps.inp"
        si_data = WaterModels.parse_file(network_path; per_unit=false)
        si_data_mn = make_multinetwork(si_data)

        # Transform the data to a per-unit system.
        pu_data_mn = deepcopy(si_data_mn)
        make_per_unit!(pu_data_mn)

        # Ensure select demand quantities are correctly transformed.
        flow_min_in_si = pu_data_mn["nw"]["1"]["demand"]["2"]["flow_min"] * pu_data_mn["base_flow"]
        @test isapprox(flow_min_in_si, si_data_mn["nw"]["1"]["demand"]["2"]["flow_min"])
        flow_max_in_si = pu_data_mn["nw"]["1"]["demand"]["2"]["flow_max"] * pu_data_mn["base_flow"]
        @test isapprox(flow_max_in_si, si_data_mn["nw"]["1"]["demand"]["2"]["flow_max"])
        flow_nominal_in_si = pu_data_mn["nw"]["1"]["demand"]["2"]["flow_nominal"] * pu_data_mn["base_flow"]
        @test isapprox(flow_nominal_in_si, si_data_mn["nw"]["1"]["demand"]["2"]["flow_nominal"])        
    end
end


@testset "make_si_units! (demand)" begin
    @testset "demand, single network" begin
        # Load the data without making a per-unit transformation.
        network_path = "../test/data/epanet/multinetwork/reservoir-hw-lps.inp"
        si_data = WaterModels.parse_file(network_path; per_unit=false)

        # Transform the data to a per-unit system.
        data = deepcopy(si_data)
        make_per_unit!(data)

        # Transform per-unit data back to SI units.
        make_si_units!(data)

        # Ensure select demand quantities are correctly transformed.
        @test isapprox(data["demand"]["2"]["flow_min"], si_data["demand"]["2"]["flow_min"])
        @test isapprox(data["demand"]["2"]["flow_max"], si_data["demand"]["2"]["flow_max"])
        @test isapprox(data["demand"]["2"]["flow_nominal"], si_data["demand"]["2"]["flow_nominal"])
        @test isapprox(data["time_series"]["demand"]["2"]["flow_min"], si_data["time_series"]["demand"]["2"]["flow_min"])
        @test isapprox(data["time_series"]["demand"]["2"]["flow_max"], si_data["time_series"]["demand"]["2"]["flow_max"])
        @test isapprox(data["time_series"]["demand"]["2"]["flow_nominal"], si_data["time_series"]["demand"]["2"]["flow_nominal"])
    end

    @testset "demand, multinetwork" begin
        # Load the data without making a per-unit transformation.
        network_path = "../test/data/epanet/multinetwork/reservoir-hw-lps.inp"
        si_data = WaterModels.parse_file(network_path; per_unit=false)
        si_data_mn = make_multinetwork(si_data)

        # Transform the data to a per-unit system.
        data_mn = deepcopy(si_data_mn)
        make_per_unit!(data_mn)

        # Transform per-unit data back to SI units.
        make_si_units!(data_mn)

        # Ensure select demand quantities are correctly transformed.
        @test isapprox(data_mn["nw"]["1"]["demand"]["2"]["flow_min"], si_data_mn["nw"]["1"]["demand"]["2"]["flow_min"])
        @test isapprox(data_mn["nw"]["1"]["demand"]["2"]["flow_max"], si_data_mn["nw"]["1"]["demand"]["2"]["flow_max"])
        @test isapprox(data_mn["nw"]["1"]["demand"]["2"]["flow_nominal"], si_data_mn["nw"]["1"]["demand"]["2"]["flow_nominal"])
    end
end
