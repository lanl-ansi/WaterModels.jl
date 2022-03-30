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
        wm = instantiate_model(network_path, NCWaterModel, build_wf)
        @test isa(wm, NCWaterModel) && isa(wm, AbstractWaterModel)
    end

    @testset "instantiate_model (with network dictionary input)" begin
        wm = instantiate_model(parse_file(network_path), NCWaterModel, build_wf)
        @test isa(wm, NCWaterModel) && isa(wm, AbstractWaterModel)
    end

    @testset "_ref_add_core!" begin
        wm = instantiate_model(parse_file(network_path), NCWaterModel, build_wf)
        WaterModels._ref_add_core!(wm.ref[:it][wm_it_sym][:nw])
        @test length(ref(wm, :pipe)) == 1
    end

    @testset "ref_add_core!" begin
        wm = instantiate_model(parse_file(network_path), NCWaterModel, build_wf)
        ref_add_core!(wm.ref)
        @test length(ref(wm, :pipe)) == 1
    end

    @testset "run_model (with file path input)" begin
        result = run_model(
            network_path,
            NCWaterModel,
            nlp_solver,
            build_wf;
            relax_integrality = true,
        )

        @test _is_valid_status(result["termination_status"])
    end

    @testset "solve_model (with non-matching multinetwork)" begin
        data, type = parse_file(network_path), NCWaterModel

        @test_throws ErrorException solve_model(
            data,
            type,
            nlp_solver,
            build_wf;
            multinetwork = true,
            relax_integrality = true,
        )
    end

    @testset "solve_model (with file path input)" begin
        result = solve_model(
            network_path,
            NCWaterModel,
            nlp_solver,
            build_wf;
            relax_integrality = true,
        )

        @test _is_valid_status(result["termination_status"])
    end

    @testset "solve_model (with network dictionary input)" begin
        result = solve_model(
            parse_file(network_path),
            NCWaterModel,
            nlp_solver,
            build_wf;
            relax_integrality = true,
        )

        @test _is_valid_status(result["termination_status"])
    end

    @testset "ismultinetwork helper function" begin
        wm = instantiate_model(network_path, NCWaterModel, build_wf)
        @test WaterModels.ismultinetwork(wm) == false
    end

    @testset "nw_ids helper function" begin
        wm = instantiate_model(network_path, NCWaterModel, build_wf)
        @test Vector{Int}(collect(nw_ids(wm))) == Vector{Int}([_IM.nw_id_default])
    end

    @testset "nws helper function" begin
        wm = instantiate_model(network_path, NCWaterModel, build_wf)
        @test Vector{Int}(collect(keys(nws(wm)))) == Vector{Int}([_IM.nw_id_default])
    end

    @testset "ids helper functions" begin
        wm = instantiate_model(network_path, NCWaterModel, build_wf)
        @test Vector{Int}(collect(ids(wm, _IM.nw_id_default, :pipe))) == Vector{Int}([1])
        @test Vector{Int}(collect(ids(wm, :pipe))) == Vector{Int}([1])
    end

    @testset "ref helper functions" begin
        wm = instantiate_model(network_path, NCWaterModel, build_wf)
        @test haskey(ref(wm, _IM.nw_id_default)[:pipe], 1)
        @test haskey(ref(wm, _IM.nw_id_default, :pipe), 1)
        @test haskey(ref(wm, _IM.nw_id_default, :pipe, 1), "diameter")
        @test isa(ref(wm, _IM.nw_id_default, :pipe, 1, "diameter"), Float64)
        @test haskey(ref(wm, :pipe), 1)
        @test haskey(ref(wm, :pipe, 1), "diameter")
        @test isa(ref(wm, :pipe, 1, "diameter"), Float64)
    end

    @testset "var helper functions" begin
        wm = instantiate_model(network_path, NCWaterModel, build_wf)
        @test JuMP.is_valid(wm.model, var(wm, _IM.nw_id_default)[:q_pipe][1])
        @test JuMP.is_valid(wm.model, var(wm, _IM.nw_id_default, :q_pipe)[1])
        @test JuMP.is_valid(wm.model, var(wm, _IM.nw_id_default, :q_pipe, 1))
        @test JuMP.is_valid(wm.model, var(wm, :q_pipe)[1])
        @test JuMP.is_valid(wm.model, var(wm, :q_pipe, 1))
    end

    @testset "con helper functions" begin
        wm = instantiate_model(network_path, NCWaterModel, build_wf)
        @test JuMP.is_valid(wm.model, con(wm, _IM.nw_id_default)[:flow_conservation][1])
        @test JuMP.is_valid(wm.model, con(wm, _IM.nw_id_default, :flow_conservation)[1])
        @test JuMP.is_valid(wm.model, con(wm, _IM.nw_id_default, :flow_conservation, 1))
        @test JuMP.is_valid(wm.model, con(wm, :flow_conservation)[1])
        @test JuMP.is_valid(wm.model, con(wm, :flow_conservation, 1))
    end

    @testset "sol helper functions" begin
        wm = instantiate_model(network_path, NCWaterModel, build_wf)
        @test haskey(sol(wm, _IM.nw_id_default)[:pipe], 1)
        @test haskey(sol(wm)[:pipe], 1)
    end
end
