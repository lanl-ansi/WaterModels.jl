function _choose_solver(wm::AbstractWaterModel, nl_solver::JuMP.MOI.OptimizerWithAttributes, mip_solver::JuMP.MOI.OptimizerWithAttributes)
    if isa(wm, AbstractNonlinearModel)
        f = Juniper.register(head_loss_args(wm)..., autodiff=false)
        return JuMP.optimizer_with_attributes(Juniper.Optimizer,
            "nl_solver"=>nl_solver, "registered_functions"=>[f], "log_levels"=>[],
            "allow_almost_solved_integral" => false)
    else
        return mip_solver
    end
end


function _build_null_model(wm::AbstractWaterModel)
    # No model variables or constraints will be added.
end


function _is_valid_status(status::JuMP.TerminationStatusCode)
    return status in [OPTIMAL, LOCALLY_SOLVED, ALMOST_LOCALLY_SOLVED]
end