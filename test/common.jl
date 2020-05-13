function make_juniper_optimizer(wm::AbstractWaterModel, nl_solver::_MOI.AbstractOptimizer)
    f = Juniper.register(head_loss_args(wm)..., autodiff=false)
    return JuMP.optimizer_with_attributes(Juniper.Optimizer,
        "nl_solver"=>nl_solver, "registered_functions"=>[f], "log_levels"=>[])
end
