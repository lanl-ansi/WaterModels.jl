@testset "Optimal Water Flow Problems" begin
    @testset "Reduced van Zyl network, multinetwork, MICP formulation." begin
        data = parse_file("../test/data/epanet/van_zyl-3_steps.inp")
        mn_data = WaterModels.make_multinetwork(data)
        ext = Dict(:pump_breakpoints=>3)
        wm = instantiate_model(mn_data, MICPWaterModel, WaterModels.build_mn_owf_agm, ext=ext)
        f = Juniper.register(head_loss_args(wm)..., autodiff=false)
        juniper = JuMP.optimizer_with_attributes(Juniper.Optimizer,
            "nl_solver"=>ipopt, "registered_functions"=>[f], "log_levels"=>[])
        ## TODO: The below takes too long to execute.
        #solution = _IM.optimize_model!(wm, optimizer=juniper)
        #@test solution["termination_status"] == LOCALLY_SOLVED
    end

    @testset "Reduced van Zyl network, multinetwork, MILP formulation." begin
        data = parse_file("../test/data/epanet/van_zyl-3_steps.inp")
        mn_data = WaterModels.make_multinetwork(data)
        ext = Dict(:pipe_breakpoints=>5, :pump_breakpoints=>5)
        wm = instantiate_model(mn_data, MILPWaterModel, WaterModels.build_mn_owf_agm, ext=ext)
        solution = _IM.optimize_model!(wm, optimizer=cbc)
        @test solution["termination_status"] == OPTIMAL
    end

    @testset "Reduced van Zyl network, multinetwork, MILP-R formulation." begin
        data = parse_file("../test/data/epanet/van_zyl-3_steps.inp")
        mn_data = WaterModels.make_multinetwork(data)
        ext = Dict(:pipe_breakpoints=>3, :pump_breakpoints=>3)
        wm = instantiate_model(mn_data, MILPRWaterModel, WaterModels.build_mn_owf_agm, ext=ext)
        solution = _IM.optimize_model!(wm, optimizer=cbc)
        @test solution["termination_status"] == OPTIMAL
    end

    @testset "Reduced van Zyl network, multinetwork, NCNLP formulation." begin
        data = parse_file("../test/data/epanet/van_zyl-3_steps.inp")
        mn_data = WaterModels.make_multinetwork(data)
        wm = instantiate_model(mn_data, NCNLPWaterModel, WaterModels.build_mn_owf_agm)
        f = Juniper.register(head_loss_args(wm)..., autodiff=false)
        juniper = JuMP.optimizer_with_attributes(Juniper.Optimizer,
            "nl_solver"=>ipopt, "registered_functions"=>[f], "log_levels"=>[])
        solution = _IM.optimize_model!(wm, optimizer=juniper)
        @test solution["termination_status"] == LOCALLY_SOLVED
    end
end
