function _make_juniper(wm::AbstractWaterModel, nl_solver::_MOI.OptimizerWithAttributes)
    f = Juniper.register(head_loss_args(wm)..., autodiff=false)
    return JuMP.optimizer_with_attributes(Juniper.Optimizer,
        "nl_solver"=>nl_solver, "registered_functions"=>[f], "log_levels"=>[])
end
