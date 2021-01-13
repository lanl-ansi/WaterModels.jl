function _make_juniper(wm::AbstractWaterModel, nl_solver::_MOI.OptimizerWithAttributes; register::Bool=true)
    if isa(wm, AbstractNonlinearModel)
        f = Juniper.register(head_loss_args(wm)..., autodiff=false)
        return JuMP.optimizer_with_attributes(Juniper.Optimizer,
            "nl_solver"=>nl_solver, "registered_functions"=>[f], "log_levels"=>[],
            "allow_almost_solved_integral" => false)
    else
        return JuMP.optimizer_with_attributes(Juniper.Optimizer,
            "nl_solver"=>nl_solver, "log_levels"=>[],
            "allow_almost_solved_integral" => false)
    end
end
