@testset "src/core/data.jl" begin
    @testset "make_per_unit!" begin
        @testset "metadata" begin
            # Load the data without making a per-unit transformation.
            si_data = WaterModels.parse_file("../examples/data/json/shamir.json"; per_unit=false)
            @test si_data["per_unit"] == false

            # Transform the data to a per-unit system.
            pu_data = deepcopy(si_data)
            make_per_unit!(pu_data)
            @test pu_data["per_unit"] == true
        end

        @testset "node, single network" begin
            # Load the data without making a per-unit transformation.
            network_path = "../test/data/epanet/multinetwork/reservoir-hw-lps.inp"
            si_data = WaterModels.parse_file(network_path; per_unit=false)

            # Transform the data to a per-unit system.
            pu_data = deepcopy(si_data)
            make_per_unit!(pu_data)

            # Ensure select node quantities are correctly transformed.
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

        @testset "tank, single network" begin
            # Load the data without making a per-unit transformation.
            network_path = "../test/data/epanet/multinetwork/tank-hw-lps.inp"
            si_data = WaterModels.parse_file(network_path; per_unit=false)
            si_data["tank"]["1"]["min_vol"] = 2.0
    
            # Transform the data to a per-unit system.
            pu_data = deepcopy(si_data)
            make_per_unit!(pu_data)
    
            # Ensure select tank quantities are correctly transformed.
            min_level_in_si = pu_data["tank"]["1"]["min_level"] * pu_data["base_length"]
            @test isapprox(min_level_in_si, si_data["tank"]["1"]["min_level"])
            max_level_in_si = pu_data["tank"]["1"]["max_level"] * pu_data["base_length"]
            @test isapprox(max_level_in_si, si_data["tank"]["1"]["max_level"])
            init_level_in_si = pu_data["tank"]["1"]["init_level"] * pu_data["base_length"]
            @test isapprox(init_level_in_si, si_data["tank"]["1"]["init_level"])
            diameter_in_si = pu_data["tank"]["1"]["diameter"] * pu_data["base_length"]
            @test isapprox(diameter_in_si, si_data["tank"]["1"]["diameter"])
            min_vol_in_si = pu_data["tank"]["1"]["min_vol"] * pu_data["base_length"]^3
            @test isapprox(min_vol_in_si, si_data["tank"]["1"]["min_vol"])
        end
    
        @testset "tank, multinetwork" begin
            # Load the data without making a per-unit transformation.
            network_path = "../test/data/epanet/multinetwork/tank-hw-lps.inp"
            si_data = WaterModels.parse_file(network_path; per_unit=false)
            si_data["tank"]["1"]["min_vol"] = 2.0
            si_data_mn = make_multinetwork(si_data)
    
            # Transform the data to a per-unit system.
            pu_data_mn = deepcopy(si_data_mn)
            make_per_unit!(pu_data_mn)
    
            # Ensure select tank quantities are correctly transformed.
            min_level_in_si = pu_data_mn["nw"]["1"]["tank"]["1"]["min_level"] * pu_data_mn["base_length"]
            @test isapprox(min_level_in_si, si_data_mn["nw"]["1"]["tank"]["1"]["min_level"])
            max_level_in_si = pu_data_mn["nw"]["1"]["tank"]["1"]["max_level"] * pu_data_mn["base_length"]
            @test isapprox(max_level_in_si, si_data_mn["nw"]["1"]["tank"]["1"]["max_level"])
            init_level_in_si = pu_data_mn["nw"]["1"]["tank"]["1"]["init_level"] * pu_data_mn["base_length"]
            @test isapprox(init_level_in_si, si_data_mn["nw"]["1"]["tank"]["1"]["init_level"])
            diameter_in_si = pu_data_mn["nw"]["1"]["tank"]["1"]["diameter"] * pu_data_mn["base_length"]
            @test isapprox(diameter_in_si, si_data_mn["nw"]["1"]["tank"]["1"]["diameter"])
            min_vol_in_si = pu_data_mn["nw"]["1"]["tank"]["1"]["min_vol"] * pu_data_mn["base_length"]^3
            @test isapprox(min_vol_in_si, si_data_mn["nw"]["1"]["tank"]["1"]["min_vol"])
        end

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

        @testset "pump, single network" begin
            # Load the data without making a per-unit transformation.
            network_path = "../test/data/epanet/snapshot/pump-hw-lps.inp"
            si_data = WaterModels.parse_file(network_path; per_unit=false)

            # Transform the data to a per-unit system.
            pu_data = deepcopy(si_data)
            make_per_unit!(pu_data)

            # Ensure select pump quantities are correctly transformed.
            flow_min_in_si = pu_data["pump"]["1"]["flow_min"] * pu_data["base_flow"]
            @test isapprox(flow_min_in_si, si_data["pump"]["1"]["flow_min"])
            flow_min_forward_in_si = pu_data["pump"]["1"]["flow_min_forward"] * pu_data["base_flow"]
            @test isapprox(flow_min_forward_in_si, si_data["pump"]["1"]["flow_min_forward"])
            flow_max_in_si = pu_data["pump"]["1"]["flow_max"] * pu_data["base_flow"]
            @test isapprox(flow_max_in_si, si_data["pump"]["1"]["flow_max"])
            flow_max_reverse_in_si = pu_data["pump"]["1"]["flow_max_reverse"] * pu_data["base_flow"]
            @test isapprox(flow_max_reverse_in_si, si_data["pump"]["1"]["flow_max_reverse"])
            head_curve_flow_in_si = [x[1] for x in pu_data["pump"]["1"]["head_curve"]] * pu_data["base_flow"]
            @test isapprox(head_curve_flow_in_si, [x[1] for x in si_data["pump"]["1"]["head_curve"]])
            head_curve_head_in_si = [x[2] for x in pu_data["pump"]["1"]["head_curve"]] * pu_data["base_head"]
            @test isapprox(head_curve_head_in_si, [x[2] for x in si_data["pump"]["1"]["head_curve"]])
        end

        @testset "pump, multinetwork" begin
            # Load the data without making a per-unit transformation.
            network_path = "../test/data/epanet/multinetwork/pump-hw-lps.inp"
            si_data = WaterModels.parse_file(network_path; per_unit=false)
            si_data_mn = make_multinetwork(si_data)

            # Transform the data to a per-unit system.
            pu_data_mn = deepcopy(si_data_mn)
            make_per_unit!(pu_data_mn)

            # Ensure select pump quantities are correctly transformed.
            flow_min_in_si = pu_data_mn["nw"]["2"]["pump"]["1"]["flow_min"] * pu_data_mn["base_flow"]
            @test isapprox(flow_min_in_si, si_data_mn["nw"]["2"]["pump"]["1"]["flow_min"])
            flow_min_forward_in_si = pu_data_mn["nw"]["2"]["pump"]["1"]["flow_min_forward"] * pu_data_mn["base_flow"]
            @test isapprox(flow_min_forward_in_si, si_data_mn["nw"]["2"]["pump"]["1"]["flow_min_forward"])
            flow_max_in_si = pu_data_mn["nw"]["2"]["pump"]["1"]["flow_max"] * pu_data_mn["base_flow"]
            @test isapprox(flow_max_in_si, si_data_mn["nw"]["2"]["pump"]["1"]["flow_max"])
            flow_max_reverse_in_si = pu_data_mn["nw"]["2"]["pump"]["1"]["flow_max_reverse"] * pu_data_mn["base_flow"]
            @test isapprox(flow_max_reverse_in_si, si_data_mn["nw"]["2"]["pump"]["1"]["flow_max_reverse"])
            head_curve_flow_in_si = [x[1] for x in pu_data_mn["nw"]["2"]["pump"]["1"]["head_curve"]] * pu_data_mn["base_flow"]
            @test isapprox(head_curve_flow_in_si, [x[1] for x in si_data_mn["nw"]["2"]["pump"]["1"]["head_curve"]])
            head_curve_head_in_si = [x[2] for x in pu_data_mn["nw"]["2"]["pump"]["1"]["head_curve"]] * pu_data_mn["base_head"]
            @test isapprox(head_curve_head_in_si, [x[2] for x in si_data_mn["nw"]["2"]["pump"]["1"]["head_curve"]])
        end

        @testset "regulator, single network" begin
            # Load the data without making a per-unit transformation.
            network_path = "../test/data/epanet/snapshot/prv-hw-lps.inp"
            si_data = WaterModels.parse_file(network_path; per_unit=false)
            si_data["regulator"]["2"]["minor_loss"] = 1.0
            si_data["regulator"]["2"]["flow_min_forward"] = 0.5

            # Transform the data to a per-unit system.
            pu_data = deepcopy(si_data)
            make_per_unit!(pu_data)

            # Ensure select regulator quantities are correctly transformed.
            diameter_in_si = pu_data["regulator"]["2"]["diameter"] * pu_data["base_length"]
            @test isapprox(diameter_in_si, si_data["regulator"]["2"]["diameter"])
            flow_min_in_si = pu_data["regulator"]["2"]["flow_min"] * pu_data["base_flow"]
            @test isapprox(flow_min_in_si, si_data["regulator"]["2"]["flow_min"])
            flow_min_forward_in_si = pu_data["regulator"]["2"]["flow_min_forward"] * pu_data["base_flow"]
            @test isapprox(flow_min_forward_in_si, si_data["regulator"]["2"]["flow_min_forward"])
            flow_max_in_si = pu_data["regulator"]["2"]["flow_max"] * pu_data["base_flow"]
            @test isapprox(flow_max_in_si, si_data["regulator"]["2"]["flow_max"])
            flow_max_reverse_in_si = pu_data["regulator"]["2"]["flow_max_reverse"] * pu_data["base_flow"]
            @test isapprox(flow_max_reverse_in_si, si_data["regulator"]["2"]["flow_max_reverse"])
            setting_in_si = pu_data["regulator"]["2"]["setting"] * pu_data["base_head"]
            @test isapprox(setting_in_si, si_data["regulator"]["2"]["setting"])
        end

        @testset "regulator, multinetwork" begin
            network_path = "../test/data/epanet/multinetwork/prv-hw-lps.inp"
            si_data = WaterModels.parse_file(network_path; per_unit=false)
            si_data["regulator"]["2"]["minor_loss"] = 1.0
            si_data["regulator"]["2"]["flow_min_forward"] = 0.5
            si_data_mn = make_multinetwork(si_data)

            # Transform the data to a per-unit system.
            pu_data_mn = deepcopy(si_data_mn)
            make_per_unit!(pu_data_mn)

            # Ensure select regulator quantities are correctly transformed.
            diameter_in_si = pu_data_mn["nw"]["2"]["regulator"]["2"]["diameter"] * pu_data_mn["base_length"]
            @test isapprox(diameter_in_si, si_data_mn["nw"]["2"]["regulator"]["2"]["diameter"])
            flow_min_in_si = pu_data_mn["nw"]["2"]["regulator"]["2"]["flow_min"] * pu_data_mn["base_flow"]
            @test isapprox(flow_min_in_si, si_data_mn["nw"]["2"]["regulator"]["2"]["flow_min"])
            flow_min_forward_in_si = pu_data_mn["nw"]["2"]["regulator"]["2"]["flow_min_forward"] * pu_data_mn["base_flow"]
            @test isapprox(flow_min_forward_in_si, si_data_mn["nw"]["2"]["regulator"]["2"]["flow_min_forward"])
            flow_max_in_si = pu_data_mn["nw"]["2"]["regulator"]["2"]["flow_max"] * pu_data_mn["base_flow"]
            @test isapprox(flow_max_in_si, si_data_mn["nw"]["2"]["regulator"]["2"]["flow_max"])
            flow_max_reverse_in_si = pu_data_mn["nw"]["2"]["regulator"]["2"]["flow_max_reverse"] * pu_data_mn["base_flow"]
            @test isapprox(flow_max_reverse_in_si, si_data_mn["nw"]["2"]["regulator"]["2"]["flow_max_reverse"])
            minor_loss_in_si = pu_data_mn["nw"]["2"]["regulator"]["2"]["minor_loss"] * pu_data_mn["base_flow"]
            @test isapprox(minor_loss_in_si, si_data_mn["nw"]["2"]["regulator"]["2"]["minor_loss"])
            setting_in_si = pu_data_mn["nw"]["2"]["regulator"]["2"]["setting"] * pu_data_mn["base_head"]
            @test isapprox(setting_in_si, si_data_mn["nw"]["2"]["regulator"]["2"]["setting"])
        end

        @testset "short pipe, single network" begin
            network_path = "../test/data/epanet/snapshot/short-pipe-lps.inp"
            si_data = WaterModels.parse_file(network_path; per_unit=false)
            WaterModels.convert_short_pipes!(si_data)
            si_data["short_pipe"]["1"]["minor_loss"] = 1.0
            si_data["short_pipe"]["1"]["flow_max_reverse"] = -0.5
            si_data["short_pipe"]["1"]["flow_min_forward"] = 0.5

            # Transform the data to a per-unit system.
            pu_data = deepcopy(si_data)
            make_per_unit!(pu_data)

            # Ensure select short pipe quantities are correctly transformed.
            flow_min_in_si = pu_data["short_pipe"]["1"]["flow_min"] * pu_data["base_flow"]
            @test isapprox(flow_min_in_si, si_data["short_pipe"]["1"]["flow_min"])
            flow_min_forward_in_si = pu_data["short_pipe"]["1"]["flow_min_forward"] * pu_data["base_flow"]
            @test isapprox(flow_min_forward_in_si, si_data["short_pipe"]["1"]["flow_min_forward"])
            flow_max_in_si = pu_data["short_pipe"]["1"]["flow_max"] * pu_data["base_flow"]
            @test isapprox(flow_max_in_si, si_data["short_pipe"]["1"]["flow_max"])
            flow_max_reverse_in_si = pu_data["short_pipe"]["1"]["flow_max_reverse"] * pu_data["base_flow"]
            @test isapprox(flow_max_reverse_in_si, si_data["short_pipe"]["1"]["flow_max_reverse"])
            minor_loss_in_si = pu_data["short_pipe"]["1"]["minor_loss"] * pu_data["base_flow"]
            @test isapprox(minor_loss_in_si, si_data["short_pipe"]["1"]["minor_loss"])
        end

        @testset "short pipe, multinetwork" begin
            network_path = "../test/data/epanet/multinetwork/short-pipe-lps.inp"
            si_data = WaterModels.parse_file(network_path; per_unit=false)
            WaterModels.convert_short_pipes!(si_data)
            si_data["short_pipe"]["1"]["minor_loss"] = 1.0
            si_data["short_pipe"]["1"]["flow_max_reverse"] = -0.5
            si_data["short_pipe"]["1"]["flow_min_forward"] = 0.5
            si_data_mn = make_multinetwork(si_data)

            # Transform the data to a per-unit system.
            pu_data_mn = deepcopy(si_data_mn)
            make_per_unit!(pu_data_mn)

            # Ensure select short pipe quantities are correctly transformed.
            flow_min_in_si = pu_data_mn["nw"]["2"]["short_pipe"]["1"]["flow_min"] * pu_data_mn["base_flow"]
            @test isapprox(flow_min_in_si, si_data_mn["nw"]["2"]["short_pipe"]["1"]["flow_min"])
            flow_min_forward_in_si = pu_data_mn["nw"]["2"]["short_pipe"]["1"]["flow_min_forward"] * pu_data_mn["base_flow"]
            @test isapprox(flow_min_forward_in_si, si_data_mn["nw"]["2"]["short_pipe"]["1"]["flow_min_forward"])
            flow_max_in_si = pu_data_mn["nw"]["2"]["short_pipe"]["1"]["flow_max"] * pu_data_mn["base_flow"]
            @test isapprox(flow_max_in_si, si_data_mn["nw"]["2"]["short_pipe"]["1"]["flow_max"])
            flow_max_reverse_in_si = pu_data_mn["nw"]["2"]["short_pipe"]["1"]["flow_max_reverse"] * pu_data_mn["base_flow"]
            @test isapprox(flow_max_reverse_in_si, si_data_mn["nw"]["2"]["short_pipe"]["1"]["flow_max_reverse"])
            minor_loss_in_si = pu_data_mn["nw"]["2"]["short_pipe"]["1"]["minor_loss"] * pu_data_mn["base_flow"]
            @test isapprox(minor_loss_in_si, si_data_mn["nw"]["2"]["short_pipe"]["1"]["minor_loss"])
        end

        @testset "valve, single network" begin
            network_path = "../test/data/epanet/snapshot/shutoff_valve-hw-lps.inp"
            si_data = WaterModels.parse_file(network_path; per_unit=false)
            si_data["valve"]["1"]["minor_loss"] = 1.0
            si_data["valve"]["1"]["flow_max_reverse"] = -0.5
            si_data["valve"]["1"]["flow_min_forward"] = 0.5

            # Transform the data to a per-unit system.
            pu_data = deepcopy(si_data)
            make_per_unit!(pu_data)

            # Ensure select valve quantities are correctly transformed.
            flow_min_in_si = pu_data["valve"]["1"]["flow_min"] * pu_data["base_flow"]
            @test isapprox(flow_min_in_si, si_data["valve"]["1"]["flow_min"])
            flow_min_forward_in_si = pu_data["valve"]["1"]["flow_min_forward"] * pu_data["base_flow"]
            @test isapprox(flow_min_forward_in_si, si_data["valve"]["1"]["flow_min_forward"])
            flow_max_in_si = pu_data["valve"]["1"]["flow_max"] * pu_data["base_flow"]
            @test isapprox(flow_max_in_si, si_data["valve"]["1"]["flow_max"])
            flow_max_reverse_in_si = pu_data["valve"]["1"]["flow_max_reverse"] * pu_data["base_flow"]
            @test isapprox(flow_max_reverse_in_si, si_data["valve"]["1"]["flow_max_reverse"])
            minor_loss_in_si = pu_data["valve"]["1"]["minor_loss"] * pu_data["base_flow"]
            @test isapprox(minor_loss_in_si, si_data["valve"]["1"]["minor_loss"])
        end

        @testset "valve, multinetwork" begin
            network_path = "../test/data/epanet/multinetwork/shutoff_valve-hw-lps.inp"
            si_data = WaterModels.parse_file(network_path; per_unit=false)
            si_data["valve"]["1"]["minor_loss"] = 1.0
            si_data["valve"]["1"]["flow_max_reverse"] = -0.5
            si_data["valve"]["1"]["flow_min_forward"] = 0.5
            si_data_mn = make_multinetwork(si_data)

            # Transform the data to a per-unit system.
            pu_data_mn = deepcopy(si_data_mn)
            make_per_unit!(pu_data_mn)

            # Ensure select valve quantities are correctly transformed.
            flow_min_in_si = pu_data_mn["nw"]["2"]["valve"]["1"]["flow_min"] * pu_data_mn["base_flow"]
            @test isapprox(flow_min_in_si, si_data_mn["nw"]["2"]["valve"]["1"]["flow_min"])
            flow_min_forward_in_si = pu_data_mn["nw"]["2"]["valve"]["1"]["flow_min_forward"] * pu_data_mn["base_flow"]
            @test isapprox(flow_min_forward_in_si, si_data_mn["nw"]["2"]["valve"]["1"]["flow_min_forward"])
            flow_max_in_si = pu_data_mn["nw"]["2"]["valve"]["1"]["flow_max"] * pu_data_mn["base_flow"]
            @test isapprox(flow_max_in_si, si_data_mn["nw"]["2"]["valve"]["1"]["flow_max"])
            flow_max_reverse_in_si = pu_data_mn["nw"]["2"]["valve"]["1"]["flow_max_reverse"] * pu_data_mn["base_flow"]
            @test isapprox(flow_max_reverse_in_si, si_data_mn["nw"]["2"]["valve"]["1"]["flow_max_reverse"])
            minor_loss_in_si = pu_data_mn["nw"]["2"]["valve"]["1"]["minor_loss"] * pu_data_mn["base_flow"]
            @test isapprox(minor_loss_in_si, si_data_mn["nw"]["2"]["valve"]["1"]["minor_loss"])
        end
    end

    @testset "make_si_units!" begin
        data = WaterModels.parse_file("../examples/data/json/shamir.json")
        make_si_units!(data)
        @test data["per_unit"] == false

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

            # Ensure select node quantities are correctly transformed.
            @test isapprox(data_mn["nw"]["1"]["node"]["1"]["head_min"], si_data_mn["nw"]["1"]["node"]["1"]["head_min"])
            @test isapprox(data_mn["nw"]["1"]["node"]["1"]["head_max"], si_data_mn["nw"]["1"]["node"]["1"]["head_max"])
            @test isapprox(data_mn["nw"]["1"]["node"]["1"]["head_nominal"], si_data_mn["nw"]["1"]["node"]["1"]["head_nominal"])
            @test isapprox(data_mn["nw"]["1"]["node"]["1"]["elevation"], si_data_mn["nw"]["1"]["node"]["1"]["elevation"])
        end

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

        @testset "tank, single network" begin
            # Load the data without making a per-unit transformation.
            network_path = "../test/data/epanet/multinetwork/tank-hw-lps.inp"
            si_data = WaterModels.parse_file(network_path; per_unit=false)
            si_data["tank"]["1"]["min_vol"] = 2.0
    
            # Transform the data to a per-unit system.
            data = deepcopy(si_data)
            make_per_unit!(data)
    
            # Transform per-unit data back to SI units.
            make_si_units!(data)
    
            # Ensure select tank quantities are correctly transformed.
            @test isapprox(data["tank"]["1"]["min_level"], si_data["tank"]["1"]["min_level"])
            @test isapprox(data["tank"]["1"]["max_level"], si_data["tank"]["1"]["max_level"])
            @test isapprox(data["tank"]["1"]["init_level"], si_data["tank"]["1"]["init_level"])
            @test isapprox(data["tank"]["1"]["diameter"], si_data["tank"]["1"]["diameter"])
            @test isapprox(data["tank"]["1"]["min_vol"], si_data["tank"]["1"]["min_vol"])
        end
    
        @testset "tank, multinetwork" begin
            # Load the data without making a per-unit transformation.
            network_path = "../test/data/epanet/multinetwork/tank-hw-lps.inp"
            si_data = WaterModels.parse_file(network_path; per_unit=false)
            si_data["tank"]["1"]["min_vol"] = 2.0
            si_data_mn = make_multinetwork(si_data)
    
            # Transform the data to a per-unit system.
            data_mn = deepcopy(si_data_mn)
            make_per_unit!(data_mn)
    
            # Transform per-unit data back to SI units.
            make_si_units!(data_mn)
    
            # Ensure select tank quantities are correctly transformed.
            @test isapprox(data_mn["nw"]["1"]["tank"]["1"]["min_level"], si_data_mn["nw"]["1"]["tank"]["1"]["min_level"])
            @test isapprox(data_mn["nw"]["1"]["tank"]["1"]["max_level"], si_data_mn["nw"]["1"]["tank"]["1"]["max_level"])
            @test isapprox(data_mn["nw"]["1"]["tank"]["1"]["init_level"], si_data_mn["nw"]["1"]["tank"]["1"]["init_level"])
            @test isapprox(data_mn["nw"]["1"]["tank"]["1"]["diameter"], si_data_mn["nw"]["1"]["tank"]["1"]["diameter"])
            @test isapprox(data_mn["nw"]["1"]["tank"]["1"]["min_vol"], si_data_mn["nw"]["1"]["tank"]["1"]["min_vol"])
        end

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

            # Ensure select design pipe quantities are correctly transformed.
            @test isapprox(data["des_pipe"]["1"]["length"], si_data["des_pipe"]["1"]["length"])
            @test isapprox(data["des_pipe"]["1"]["diameter"], si_data["des_pipe"]["1"]["diameter"])
            @test isapprox(data["des_pipe"]["1"]["flow_min"], si_data["des_pipe"]["1"]["flow_min"])
            @test isapprox(data["des_pipe"]["1"]["flow_max"], si_data["des_pipe"]["1"]["flow_max"])
            @test isapprox(data["des_pipe"]["1"]["minor_loss"], si_data["des_pipe"]["1"]["minor_loss"])
        end

        @testset "regulator, single network" begin
            # Load the data without making a per-unit transformation.
            network_path = "../test/data/epanet/snapshot/prv-hw-lps.inp"
            si_data = WaterModels.parse_file(network_path; per_unit=false)
            si_data["regulator"]["2"]["minor_loss"] = 1.0
            si_data["regulator"]["2"]["flow_min_forward"] = 0.5
            
            # Transform the data to a per-unit system.
            data = deepcopy(si_data)
            make_per_unit!(data)

            # Transform per-unit data back to SI units.
            make_si_units!(data)

            # Ensure select regulator quantities are correctly transformed.
            @test isapprox(data["regulator"]["2"]["diameter"], si_data["regulator"]["2"]["diameter"])
            @test isapprox(data["regulator"]["2"]["flow_min"], si_data["regulator"]["2"]["flow_min"])
            @test isapprox(data["regulator"]["2"]["flow_max"], si_data["regulator"]["2"]["flow_max"])
            @test isapprox(data["regulator"]["2"]["flow_min_forward"], si_data["regulator"]["2"]["flow_min_forward"])
            @test isapprox(data["regulator"]["2"]["flow_max_reverse"], si_data["regulator"]["2"]["flow_max_reverse"])
            @test isapprox(data["regulator"]["2"]["minor_loss"], si_data["regulator"]["2"]["minor_loss"])
            @test isapprox(data["regulator"]["2"]["setting"], si_data["regulator"]["2"]["setting"])
        end

        @testset "regulator, multinetwork" begin
            # Load the data without making a per-unit transformation.
            network_path = "../test/data/epanet/multinetwork/prv-hw-lps.inp"
            si_data = WaterModels.parse_file(network_path; per_unit=false)
            si_data["regulator"]["2"]["minor_loss"] = 1.0
            si_data["regulator"]["2"]["flow_min_forward"] = 0.5
            si_data_mn = make_multinetwork(si_data)

            # Transform the data to a per-unit system.
            data_mn = deepcopy(si_data_mn)
            make_per_unit!(data_mn)

            # Transform per-unit data back to SI units.
            make_si_units!(data_mn)

            # Ensure select regulator quantities are correctly transformed.
            @test isapprox(data_mn["nw"]["2"]["regulator"]["2"]["diameter"], si_data_mn["nw"]["2"]["regulator"]["2"]["diameter"])
            @test isapprox(data_mn["nw"]["2"]["regulator"]["2"]["flow_min"], si_data_mn["nw"]["2"]["regulator"]["2"]["flow_min"])
            @test isapprox(data_mn["nw"]["2"]["regulator"]["2"]["flow_max"], si_data_mn["nw"]["2"]["regulator"]["2"]["flow_max"])
            @test isapprox(data_mn["nw"]["2"]["regulator"]["2"]["flow_min_forward"], si_data_mn["nw"]["2"]["regulator"]["2"]["flow_min_forward"])
            @test isapprox(data_mn["nw"]["2"]["regulator"]["2"]["flow_max_reverse"], si_data_mn["nw"]["2"]["regulator"]["2"]["flow_max_reverse"])
            @test isapprox(data_mn["nw"]["2"]["regulator"]["2"]["minor_loss"], si_data_mn["nw"]["2"]["regulator"]["2"]["minor_loss"])
            @test isapprox(data_mn["nw"]["2"]["regulator"]["2"]["setting"], si_data_mn["nw"]["2"]["regulator"]["2"]["setting"])
        end

        @testset "short pipe, single network" begin
            # Load the data without making a per-unit transformation.
            network_path = "../test/data/epanet/snapshot/short-pipe-lps.inp"
            si_data = WaterModels.parse_file(network_path; per_unit=false)
            WaterModels.convert_short_pipes!(si_data)
            si_data["short_pipe"]["1"]["minor_loss"] = 1.0
            si_data["short_pipe"]["1"]["flow_max_reverse"] = -0.5
            si_data["short_pipe"]["1"]["flow_min_forward"] = 0.5
            
            # Transform the data to a per-unit system.
            data = deepcopy(si_data)
            make_per_unit!(data)

            # Transform per-unit data back to SI units.
            make_si_units!(data)

            # Ensure select short pipe quantities are correctly transformed.
            @test isapprox(data["short_pipe"]["1"]["flow_min"], si_data["short_pipe"]["1"]["flow_min"])
            @test isapprox(data["short_pipe"]["1"]["flow_max"], si_data["short_pipe"]["1"]["flow_max"])
            @test isapprox(data["short_pipe"]["1"]["flow_min_forward"], si_data["short_pipe"]["1"]["flow_min_forward"])
            @test isapprox(data["short_pipe"]["1"]["flow_max_reverse"], si_data["short_pipe"]["1"]["flow_max_reverse"])
            @test isapprox(data["short_pipe"]["1"]["minor_loss"], si_data["short_pipe"]["1"]["minor_loss"])
        end

        @testset "short pipe, multinetwork" begin
            # Load the data without making a per-unit transformation.
            network_path = "../test/data/epanet/multinetwork/short-pipe-lps.inp"
            si_data = WaterModels.parse_file(network_path; per_unit=false)
            WaterModels.convert_short_pipes!(si_data)
            si_data["short_pipe"]["1"]["minor_loss"] = 1.0
            si_data["short_pipe"]["1"]["flow_max_reverse"] = -0.5
            si_data["short_pipe"]["1"]["flow_min_forward"] = 0.5
            si_data_mn = make_multinetwork(si_data)

            # Transform the data to a per-unit system.
            data_mn = deepcopy(si_data_mn)
            make_per_unit!(data_mn)

            # Transform per-unit data back to SI units.
            make_si_units!(data_mn)

            # Ensure select short pipe quantities are correctly transformed.
            @test isapprox(data_mn["nw"]["2"]["short_pipe"]["1"]["flow_min"], si_data_mn["nw"]["2"]["short_pipe"]["1"]["flow_min"])
            @test isapprox(data_mn["nw"]["2"]["short_pipe"]["1"]["flow_max"], si_data_mn["nw"]["2"]["short_pipe"]["1"]["flow_max"])
            @test isapprox(data_mn["nw"]["2"]["short_pipe"]["1"]["flow_min_forward"], si_data_mn["nw"]["2"]["short_pipe"]["1"]["flow_min_forward"])
            @test isapprox(data_mn["nw"]["2"]["short_pipe"]["1"]["flow_max_reverse"], si_data_mn["nw"]["2"]["short_pipe"]["1"]["flow_max_reverse"])
            @test isapprox(data_mn["nw"]["2"]["short_pipe"]["1"]["minor_loss"], si_data_mn["nw"]["2"]["short_pipe"]["1"]["minor_loss"])
        end

        @testset "pump, single network" begin
            # Load the data without making a per-unit transformation.
            network_path = "../test/data/epanet/snapshot/pump-hw-lps.inp"
            si_data = WaterModels.parse_file(network_path; per_unit=false)

            # Transform the data to a per-unit system.
            data = deepcopy(si_data)
            make_per_unit!(data)

            # Transform per-unit data back to SI units.
            make_si_units!(data)

            # Ensure select pump quantities are correctly transformed.
            @test isapprox(data["pump"]["1"]["flow_min"], si_data["pump"]["1"]["flow_min"])
            @test isapprox(data["pump"]["1"]["flow_min_forward"], si_data["pump"]["1"]["flow_min_forward"])
            @test isapprox(data["pump"]["1"]["flow_max"], si_data["pump"]["1"]["flow_max"])
            @test isapprox(data["pump"]["1"]["flow_max_reverse"], si_data["pump"]["1"]["flow_max_reverse"])
            head_curve_flow = [x[1] for x in data["pump"]["1"]["head_curve"]]
            @test isapprox(head_curve_flow, [x[1] for x in si_data["pump"]["1"]["head_curve"]])
            head_curve_head = [x[2] for x in data["pump"]["1"]["head_curve"]]
            @test isapprox(head_curve_head, [x[2] for x in si_data["pump"]["1"]["head_curve"]])
        end

        @testset "pump, multinetwork" begin
            # Load the data without making a per-unit transformation.
            network_path = "../test/data/epanet/multinetwork/pump-hw-lps.inp"
            si_data = WaterModels.parse_file(network_path; per_unit=false)
            si_data_mn = make_multinetwork(si_data)

            # Transform the data to a per-unit system.
            data_mn = deepcopy(si_data_mn)
            make_per_unit!(data_mn)

            # Transform per-unit data back to SI units.
            make_si_units!(data_mn)

            # Ensure select pump quantities are correctly transformed.
            @test isapprox(data_mn["nw"]["2"]["pump"]["1"]["flow_min"], si_data_mn["nw"]["2"]["pump"]["1"]["flow_min"])
            @test isapprox(data_mn["nw"]["2"]["pump"]["1"]["flow_min_forward"], si_data_mn["nw"]["2"]["pump"]["1"]["flow_min_forward"])
            @test isapprox(data_mn["nw"]["2"]["pump"]["1"]["flow_max"], si_data_mn["nw"]["2"]["pump"]["1"]["flow_max"])
            @test isapprox(data_mn["nw"]["2"]["pump"]["1"]["flow_max_reverse"], si_data_mn["nw"]["2"]["pump"]["1"]["flow_max_reverse"])
            head_curve_flow = [x[1] for x in data_mn["nw"]["2"]["pump"]["1"]["head_curve"]]
            @test isapprox(head_curve_flow, [x[1] for x in si_data_mn["nw"]["2"]["pump"]["1"]["head_curve"]])
            head_curve_head = [x[2] for x in data_mn["nw"]["2"]["pump"]["1"]["head_curve"]]
            @test isapprox(head_curve_head, [x[2] for x in si_data_mn["nw"]["2"]["pump"]["1"]["head_curve"]])
        end

        @testset "valve, single network" begin
            # Load the data without making a per-unit transformation.
            network_path = "../test/data/epanet/snapshot/shutoff_valve-hw-lps.inp"
            si_data = WaterModels.parse_file(network_path; per_unit=false)
            si_data["valve"]["1"]["minor_loss"] = 1.0
            si_data["valve"]["1"]["flow_max_reverse"] = -0.5
            si_data["valve"]["1"]["flow_min_forward"] = 0.5
            
            # Transform the data to a per-unit system.
            data = deepcopy(si_data)
            make_per_unit!(data)

            # Transform per-unit data back to SI units.
            make_si_units!(data)

            # Ensure select valve quantities are correctly transformed.
            @test isapprox(data["valve"]["1"]["flow_min"], si_data["valve"]["1"]["flow_min"])
            @test isapprox(data["valve"]["1"]["flow_max"], si_data["valve"]["1"]["flow_max"])
            @test isapprox(data["valve"]["1"]["flow_min_forward"], si_data["valve"]["1"]["flow_min_forward"])
            @test isapprox(data["valve"]["1"]["flow_max_reverse"], si_data["valve"]["1"]["flow_max_reverse"])
            @test isapprox(data["valve"]["1"]["minor_loss"], si_data["valve"]["1"]["minor_loss"])
        end

        @testset "valve, multinetwork" begin
            # Load the data without making a per-unit transformation.
            network_path = "../test/data/epanet/multinetwork/shutoff_valve-hw-lps.inp"
            si_data = WaterModels.parse_file(network_path; per_unit=false)
            si_data["valve"]["1"]["minor_loss"] = 1.0
            si_data["valve"]["1"]["flow_max_reverse"] = -0.5
            si_data["valve"]["1"]["flow_min_forward"] = 0.5
            si_data_mn = make_multinetwork(si_data)

            # Transform the data to a per-unit system.
            data_mn = deepcopy(si_data_mn)
            make_per_unit!(data_mn)

            # Transform per-unit data back to SI units.
            make_si_units!(data_mn)

            # Ensure select valve quantities are correctly transformed.
            @test isapprox(data_mn["nw"]["2"]["valve"]["1"]["flow_min"], si_data_mn["nw"]["2"]["valve"]["1"]["flow_min"])
            @test isapprox(data_mn["nw"]["2"]["valve"]["1"]["flow_max"], si_data_mn["nw"]["2"]["valve"]["1"]["flow_max"])
            @test isapprox(data_mn["nw"]["2"]["valve"]["1"]["flow_min_forward"], si_data_mn["nw"]["2"]["valve"]["1"]["flow_min_forward"])
            @test isapprox(data_mn["nw"]["2"]["valve"]["1"]["flow_max_reverse"], si_data_mn["nw"]["2"]["valve"]["1"]["flow_max_reverse"])
            @test isapprox(data_mn["nw"]["2"]["valve"]["1"]["minor_loss"], si_data_mn["nw"]["2"]["valve"]["1"]["minor_loss"])
        end
    end

    @testset "make_temporally_aggregated_multinetwork" begin
        network = WaterModels.parse_file("../test/data/epanet/multinetwork/owf-hw-lps.inp")
        network_mn = WaterModels.make_multinetwork(network)
        network_agg = make_temporally_aggregated_multinetwork(network_mn, [["1", "2"], ["3"]])

        @test isapprox(network_mn["duration"], network_agg["duration"])
        demand_old = network_mn["nw"]["1"]["demand"]["2"]["flow_nominal"] 
        @test demand_old < network_agg["nw"]["1"]["demand"]["2"]["flow_nominal"]

        network = WaterModels.parse_file("../test/data/epanet/multinetwork/prv-hw-lps.inp")
        network_mn = WaterModels.make_multinetwork(network)
        network_agg = make_temporally_aggregated_multinetwork(network_mn, [["1", "2"], ["3"]])

        @test network_mn["duration"] == network_agg["duration"]
        @test length(network_agg["nw"]) == 2
    end

    @testset "set_flow_partitions_num!" begin
        data = WaterModels.parse_file("../examples/data/epanet/van_zyl.inp")
        WaterModels.set_flow_partitions_num!(data, 5)
        @test length(data["pipe"]["1"]["flow_partition"]) == 5
        @test length(data["pump"]["1"]["flow_partition"]) == 5
    end

    @testset "set_start! (single network)" begin
        data = WaterModels.parse_file("../test/data/epanet/snapshot/pipe-hw-lps.inp")
        data["pipe"]["1"]["q"] = 1.0 # Set the flow along the pipe.
        WaterModels.set_start!(data, "pipe", "q", "q_pipe_start")
    end

    @testset "set_start! (multinetwork)" begin
        data = WaterModels.parse_file("../test/data/epanet/snapshot/pipe-hw-lps.inp")
        data["pipe"]["1"]["q"] = 1.0 # Set the flow along the pipe.
        mn_data = WaterModels.replicate(data, 3)
        WaterModels.set_start!(mn_data, "pipe", "q", "q_pipe_start")
    end

    @testset "set_warm_start! (multinetwork)" begin
        network = WaterModels.parse_file("../examples/data/epanet/van_zyl.inp")
        network_mn = WaterModels.make_multinetwork(network)
        WaterModels.set_warm_start!(network_mn)
        @test haskey(network_mn["nw"]["1"]["pipe"]["1"], "q_start")
    end

    @testset "replicate" begin
        data = WaterModels.parse_file("../test/data/epanet/snapshot/pipe-hw-lps.inp")
        mn_data = WaterModels.replicate(data, 3)
        @test length(mn_data["nw"]) == 3
    end

    @testset "make_multinetwork shamir" begin
        network_data = WaterModels.parse_file("../test/data/epanet/multinetwork/shamir-ts.inp")

        @test !_IM.ismultinetwork(network_data)
        @test haskey(network_data, "time_series")
        base_flow = WaterModels._calc_flow_per_unit_transform(network_data)(0.02777)
        @test isapprox(network_data["time_series"]["demand"]["2"]["flow_nominal"][1], base_flow, rtol=1.0e-4)
        @test isapprox(network_data["time_series"]["demand"]["2"]["flow_nominal"][2], 0.5*base_flow, rtol=1.0e-4)
        @test isapprox(network_data["time_series"]["demand"]["2"]["flow_nominal"][3], 0.25*base_flow, rtol=1.0e-4)

        ts_length = network_data["time_series"]["num_steps"]
        mn_data = WaterModels.make_multinetwork(network_data)

        @test _IM.ismultinetwork(mn_data)
        @test !haskey(mn_data, "time_series")
        @test length(mn_data["nw"]) == ts_length
        @test isapprox(mn_data["nw"]["1"]["demand"]["2"]["flow_nominal"], base_flow, rtol=1.0e-4)
        @test isapprox(mn_data["nw"]["2"]["demand"]["2"]["flow_nominal"], 0.5*base_flow, rtol=1.0e-4)
        @test isapprox(mn_data["nw"]["3"]["demand"]["2"]["flow_nominal"], 0.25*base_flow, rtol=1.0e-4)
    end

    @testset "split_multinetwork" begin
        data = WaterModels.parse_file("../test/data/epanet/snapshot/pipe-hw-lps.inp")
        mn_data = WaterModels.replicate(data, 3)
        new_mns = WaterModels.split_multinetwork(mn_data, [["1", "2"], ["2", "3"]])

        @test sort(collect(keys(new_mns[1]["nw"]))) == ["1", "2"]
        @test sort(collect(keys(new_mns[2]["nw"]))) == ["2", "3"]
    end

    @testset "epanet_to_watermodels!" begin
        network_data = WaterModels.parse_epanet("../test/data/epanet/snapshot/short-pipe-lps.inp")
        WaterModels.epanet_to_watermodels!(network_data)
        WaterModels.correct_network_data!(network_data)
        WaterModels.convert_short_pipes!(network_data)
        @test haskey(network_data["short_pipe"], "1")

        network_data = WaterModels.parse_epanet("../test/data/epanet/snapshot/short-pipe-valve-lps.inp")
        WaterModels.epanet_to_watermodels!(network_data)
        WaterModels.correct_network_data!(network_data)
        @test haskey(network_data["valve"], "1")
    end

    @testset "_read_status!" begin
        network_data = WaterModels.parse_epanet("../test/data/epanet/snapshot/shutoff_valve-status-hw-lps.inp")
        @test network_data["pipe"]["1"]["has_valve"]
        WaterModels.epanet_to_watermodels!(network_data)
        @test haskey(network_data["valve"], "1")
    end

    @testset "_read_controls!" begin
        network_data = WaterModels.parse_epanet("../test/data/epanet/snapshot/shutoff_valve-controls-hw-lps.inp")
        @test network_data["pipe"]["1"]["has_valve"]
        WaterModels.epanet_to_watermodels!(network_data)
        @test haskey(network_data["valve"], "1")
    end

    @testset "_remove_last_networks! (multinetwork)" begin
        data = WaterModels.parse_file("../test/data/epanet/snapshot/pipe-hw-lps.inp")
        mn_data = WaterModels.replicate(data, 3)
        WaterModels._remove_last_networks!(mn_data; last_num_steps = 1)
        @test length(mn_data["nw"]) == 2
    end
end
