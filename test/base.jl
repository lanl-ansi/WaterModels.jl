@testset "src/core/base.jl" begin
    network_path = "../test/data/epanet/snapshot/pipe-hw-lps.inp"

    @testset "silence" begin
        # This should silence everything except error messages.
        WaterModels.silence()

        wm_logger = Memento.getlogger(_IM)
        @test Memento.getlevel(wm_logger) == "error"
        Memento.warn(wm_logger, "Silenced message should not be displayed.")

        wm_logger = Memento.getlogger(WaterModels)
        @test Memento.getlevel(wm_logger) == "error"
        Memento.warn(wm_logger, "Silenced message should not be displayed.")
    end

    @testset "build_ref" begin
        ref = build_ref(parse_file(network_path))
        @test haskey(ref[:it][wm_it_sym], :head_loss)
        @test haskey(ref[:it][wm_it_sym][:nw][0][:pipe], 1)
    end

    @testset "instantiate_model (with file path input)" begin
        wm = instantiate_model(network_path, LAWaterModel, build_wf)
        @test isa(wm, LAWaterModel) && isa(wm, AbstractWaterModel)
    end

    @testset "instantiate_model (with network dictionary input)" begin
        wm = instantiate_model(parse_file(network_path), LAWaterModel, build_wf)
        @test isa(wm, LAWaterModel) && isa(wm, AbstractWaterModel)
    end

    @testset "_ref_add_core!" begin
        wm = instantiate_model(parse_file(network_path), LAWaterModel, build_wf)
        WaterModels._ref_add_core!(wm.ref[:it][wm_it_sym][:nw], wm.ref[:it][wm_it_sym][:head_loss])
        @test length(ref(wm, :pipe)) == 1
    end

    @testset "ref_add_core!" begin
        wm = instantiate_model(parse_file(network_path), LAWaterModel, build_wf)
        ref_add_core!(wm.ref)
        @test length(ref(wm, :pipe)) == 1
    end

    @testset "solve_model (with non-matching multinetwork)" begin
        data, type = parse_file(network_path), LAWaterModel
        @test_throws ErrorException solve_model(data, type, cbc, build_wf; multinetwork=true)
    end

    @testset "solve_model (with file path input)" begin
        result = solve_model(network_path, LAWaterModel, cbc, build_wf)
        @test result["termination_status"] == OPTIMAL
    end

    @testset "solve_model (with network dictionary input)" begin
        result = solve_model(parse_file(network_path), LAWaterModel, cbc, build_wf)
        @test result["termination_status"] == OPTIMAL
    end
end
