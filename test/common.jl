function _choose_solver(wm::AbstractWaterModel, nl_solver::_MOI.OptimizerWithAttributes, mip_solver::_MOI.OptimizerWithAttributes)
    if isa(wm, AbstractNonlinearModel)
        f = Juniper.register(head_loss_args(wm)..., autodiff=false)
        return JuMP.optimizer_with_attributes(Juniper.Optimizer,
            "nl_solver"=>nl_solver, "registered_functions"=>[f], "log_levels"=>[],
            "allow_almost_solved_integral" => false)
    else
        return mip_solver
    end
end

function _is_valid_status(status::_MOI.TerminationStatusCode)
    return status in [OPTIMAL, LOCALLY_SOLVED, ALMOST_LOCALLY_SOLVED]
end