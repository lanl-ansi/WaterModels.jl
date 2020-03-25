#""
#function build_result(wm::AbstractWaterModel, solve_time; solution_processors=[])
#    # TODO replace with JuMP.result_count(wm.model) after version v0.21
#    # try-catch is needed until solvers reliably support ResultCount()
#    result_count = 1
#
#    try
#        result_count = _MOI.get(wm.model, _MOI.ResultCount())
#    catch
#        Memento.warn(_LOGGER, "the given optimizer does not provide the "
#            * "ResultCount() attribute, assuming the solver returned a "
#            * "solution, which may be incorrect.")
#    end
#
#    sol = Dict{String,Any}()
#
#    if result_count > 0
#        sol = build_solution(wm, post_processors=solution_processors)
#    else
#        Memento.warn(_LOGGER, "model has no results, solution cannot be built.")
#    end
#
#    data = Dict{String,Any}("name" => wm.data["name"])
#
#    if InfrastructureModels.ismultinetwork(wm.data)
#        data_nws = data["nw"] = Dict{String,Any}()
#
#        for (n, nw_data) in wm.data["nw"]
#            link_count = length(nw_data["pipe"]) + length(nw_data["pump"])
#            node_count = length(nw_data["node"])
#
#            data_nws[n] = Dict(
#                "name" => get(nw_data, "name", "anonymous"),
#                "link_count" => link_count,
#                "node_count" => node_count)
#        end
#    else
#        data["link_count"] = length(wm.data["pipe"]) + length(wm.data["pump"])
#        data["node_count"] = length(wm.data["node"])
#    end
#
#    solution = Dict{String,Any}(
#        "optimizer" => JuMP.solver_name(wm.model),
#        "termination_status" => JuMP.termination_status(wm.model),
#        "primal_status" => JuMP.primal_status(wm.model),
#        "dual_status" => JuMP.dual_status(wm.model),
#        "objective" => _guard_objective_value(wm.model),
#        "objective_lb" => _guard_objective_bound(wm.model),
#        "solve_time" => solve_time,
#        "solution" => sol,
#        "machine" => Dict("cpu" => Sys.cpu_info()[1].model,
#            "memory" => string(Sys.total_memory()/2^30, " Gb")),
#        "data" => data)
#
#    return solution
#end
#
#
#""
#function _guard_objective_value(model)
#    obj_val = NaN
#
#    try
#        obj_val = JuMP.objective_value(model)
#    catch
#    end
#
#    return obj_val
#end
#
#
#""
#function _guard_objective_bound(model)
#    obj_lb = -Inf
#
#    try
#        obj_lb = JuMP.objective_bound(model)
#    catch
#    end
#
#    return obj_lb
#end
#
#""
#function build_solution(wm::AbstractWaterModel; post_processors=[])
#    # TODO @assert that the model is solved
#    sol = _build_solution_values(wm.sol)
#    sol["per_unit"] = wm.data["per_unit"]
#
#    if ismultinetwork(wm)
#        sol["multinetwork"] = true
#    else
#        for (k,v) in sol["nw"]["$(wm.cnw)"]
#            sol[k] = v
#        end
#
#        delete!(sol, "nw")
#    end
#
#    for post_processor in post_processors
#        post_processor(wm, sol)
#    end
#
#    return sol
#end
#
#""
#function _build_solution_values(var::Dict)
#    sol = Dict{String,Any}()
#
#    for (key, val) in var
#        sol[string(key)] = _build_solution_values(val)
#    end
#
#    return sol
#end
#
#""
#function _build_solution_values(var::Array{<:Any,1})
#    return [_build_solution_values(val) for val in var]
#end
#
#""
#function _build_solution_values(var::Array{<:Any,2})
#    return [_build_solution_values(var[i,j]) for i in 1:size(var,1), j in 1:size(var,2)]
#end
#
#"Solution building support for symmetric JuMP matrix variables."
#function _build_solution_values(var::LinearAlgebra.Symmetric{T,Array{T,2}}) where T
#    return [_build_solution_values(var[i,j]) for i in 1:size(var,1), j in 1:size(var,2)]
#end
#
#""
#function _build_solution_values(var::Number)
#    return var
#end
#
#""
#function _build_solution_values(var::JuMP.VariableRef)
#    return JuMP.value(var)
#end
#
#""
#function _build_solution_values(var::JuMP.GenericAffExpr)
#    return JuMP.value(var)
#end
#
#""
#function _build_solution_values(var::JuMP.GenericQuadExpr)
#    return JuMP.value(var)
#end
#
#""
#function _build_solution_values(var::JuMP.NonlinearExpression)
#    return JuMP.value(var)
#end
#
#""
#function _build_solution_values(var::JuMP.ConstraintRef)
#    return -JuMP.dual(var)
#end
#
#""
#function _build_solution_values(var::Any)
#    Memento.warn(_LOGGER, "_build_solution_values found unknown type $(typeof(var))")
#    return var
#end
#
#"Converts the solution data into the data model's standard space, polar voltages and rectangular power"
#function sol_data_model!(wm::AbstractWaterModel, solution::Dict)
#    Memento.warn(_LOGGER, "sol_data_model! not defined for power model of type $(typeof(wm))")
#end
