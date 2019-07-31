@testset "Optimal Water Flow Problems" begin
    van_zyl_path = "../test/data/epanet/van_zyl.inp"

    @testset "van Zyl network, multinetwork, NCNLP formulation." begin
        network_data = WaterModels.parse_file(van_zyl_path)
        mn_data = WaterModels.make_multinetwork(network_data)
        wm = build_generic_model(mn_data, NCNLPWaterModel, WaterModels.post_mn_owf, multinetwork=true)
        f = Juniper.register(fun(wm, :head_loss)..., autodiff=false)
        juniper = JuMP.with_optimizer(Juniper.Optimizer, nl_solver=ipopt, registered_functions=[f], log_levels=[])
        solution = WaterModels.solve_generic_model(wm, juniper)
    end
end
