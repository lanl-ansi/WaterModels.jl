@testset "Optimal Water Flow Problems" begin
    @testset "van Zyl network, multinetwork, NCNLP formulation." begin
        data = parse_file("../test/data/epanet/van_zyl-3_steps.inp")
        mn_data = WaterModels.make_multinetwork(data)
        wm = instantiate_model(mn_data, NCNLPWaterModel, WaterModels.build_mn_owf)
        f = Juniper.register(head_loss_args(wm)..., autodiff=false)
        juniper = JuMP.optimizer_with_attributes(Juniper.Optimizer,
            "nl_solver"=>ipopt, "registered_functions"=>[f], "log_levels"=>[])
        solution = _IM.optimize_model!(wm, optimizer=juniper)
        @test solution["termination_status"] == LOCALLY_INFEASIBLE
    end
end
