@testset "Optimal Water Flow Problems" begin
    @testset "van Zyl network, multinetwork, NCNLP formulation." begin
        network_data = WaterModels.parse_file("../test/data/epanet/van_zyl.inp")
        mn_data = WaterModels.make_multinetwork(network_data)
        wm = build_model(mn_data, NCNLPWaterModel, WaterModels.post_mn_owf, multinetwork=true)
        f = Juniper.register(fun(wm, :head_loss)..., autodiff=false)
        juniper = JuMP.optimizer_with_attributes(Juniper.Optimizer, "nl_solver"=>ipopt, "registered_functions"=>[f], "log_levels"=>[])
        solution = WaterModels.optimize_model!(wm, juniper) # Currently infeasible.
    end
end
