function _choose_solver(wm::AbstractWaterModel, nl_solver::JuMP.MOI.OptimizerWithAttributes, mip_solver::JuMP.MOI.OptimizerWithAttributes)
    if isa(wm, AbstractNonlinearModel)
        return JuMP.optimizer_with_attributes(Juniper.Optimizer,
            "nl_solver"=>nl_solver, "log_levels"=>[],
            "allow_almost_solved_integral" => false)
    else
        return mip_solver
    end
end

function _is_valid_status(status::JuMP.TerminationStatusCode)
    return status in [OPTIMAL, LOCALLY_SOLVED, ALMOST_LOCALLY_SOLVED]
end
