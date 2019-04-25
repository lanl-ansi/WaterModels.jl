""
function build_solution(wm::GenericWaterModel, status, solve_time; objective = NaN, solution_builder = get_solution)
    if status != :Error
        objective = JuMP.objective_value(wm.model)
        status = JuMP.primal_status(wm.model)
    end

    sol = init_solution(wm)
    data = Dict{String,Any}("title" => wm.data["title"])

    if InfrastructureModels.ismultinetwork(wm.data)
        sol["multinetwork"] = true
        sol_nws = sol["nw"] = Dict{String,Any}()
        data_nws = data["nw"] = Dict{String,Any}()

        for (n, nw_data) in wm.data["nw"]
            sol_nw = sol_nws[n] = Dict{String,Any}()
            wm.cnw = parse(Int, n)
            solution_builder(wm, sol_nw)
            data_nws[n] = Dict("name" => get(nw_data, "name", "anonymous"),
                               "link_count" => length(wm[:ref][n][:links]),
                               "node_count" => length(wm[:ref][n][:nodes]))
        end
    else
        solution_builder(wm, sol)
        data["link_count"] = length(wm.ref[:nw][wm.cnw][:links])
        data["node_count"] = length(wm.ref[:nw][wm.cnw][:nodes])
    end

    cpu = Sys.cpu_info()[1].model
    memory = Sys.total_memory()

    solution = Dict{String,Any}(
        "termination_status" => JuMP.termination_status(wm.model),
        "primal_status" => JuMP.primal_status(wm.model),
        "dual_status" => JuMP.dual_status(wm.model),
        "objective_value" => JuMP.objective_value(wm.model),
        "solve_time" => solve_time,
        "solution" => sol,
        "machine" => Dict("cpu" => cpu, "memory" => memory),
        "data" => data)

    wm.solution = solution
    println(solution)
    return solution
end

""
function init_solution(wm::GenericWaterModel)
    data_keys = ["per_unit"]
    return Dict{String,Any}(key => wm.data[key] for key in data_keys)
end

""
function get_solution(wm::GenericWaterModel, sol::Dict{String,<:Any})
    add_pipe_flow_setpoint(sol, wm)
    #add_bus_voltage_setpoint(sol, wm)
    #add_generator_power_setpoint(sol, wm)
    #add_storage_setpoint(sol, wm)
    #add_branch_flow_setpoint(sol, wm)
    #add_dcline_flow_setpoint(sol, wm)
    #add_kcl_duals(sol, wm)
    #add_sm_duals(sol, wm) # Adds the duals of the transmission lines' thermal limits.
end

function add_pipe_flow_setpoint(sol, wm::GenericWaterModel)
    add_setpoint(sol, wm, "pipes", "q", :q)
end

"adds values based on JuMP variables"
function add_setpoint(sol::Dict{String, Any}, wm::GenericWaterModel,
                      dict_name::String, param_name::String,
                      variable_symbol::Symbol; index_name::String="id",
                      default_value=(item) -> NaN, scale=(x, item) -> x,
                      extract_var=(var, idx, item) -> var[idx],
                      sol_dict=get(sol, dict_name, Dict{String, Any}()))
    if InfrastructureModels.ismultinetwork(wm.data)
        data_dict = wm.data["nw"]["$(wm.cnw)"][dict_name]
    else
        data_dict = wm.data[dict_name]
    end

    if length(data_dict) > 0
        sol[dict_name] = sol_dict
    end

    for (i, item) in data_dict
        idx = parse(Int, item[index_name])
        sol_item = sol_dict[i] = get(sol_dict, i, Dict{String, Any}())
        sol_item[param_name] = default_value(item)

        try
            variable = extract_var(var(wm, wm.cnw, variable_symbol), idx, item)
            sol_item[param_name] = scale(JuMP.value(variable), item, 1)
        catch
        end

        sol_item[param_name] = sol_item[param_name][1]
    end
end

#function build_solution(wm::GenericWaterModel, status, solve_time;
#                        objective = NaN, solution_builder = get_solution)
#    if status != :Error
#        objective = getobjectivevalue(wm.model)
#        status = solver_status_dict(Symbol(typeof(wm.model.solver).name.module), status)
#    end
#
#    sol = init_solution(wm)
#    data = Dict{String, Any}("name" => wm.data["title"])
#
#    if wm.data["multinetwork"]
#        sol_nws = sol["nw"] = Dict{String, Any}()
#        data_nws = data["nw"] = Dict{String, Any}()
#
#        for (n, nw_data) in wm.data["nw"]
#            sol_nw = sol_nws[n] = Dict{String, Any}()
#            wm.cnw = parse(Int, n)
#            data_nws[n] = Dict(
#                "name" => get(nw_data, "name", "anonymous"),
#                "node_count" => length(nw_data["junctions"]) + length(nw_data["reservoirs"]),
#                "link_count" => length(nw_data["pipes"])
#            )
#
#            if status != :LocalInfeasible
#                solution_builder(wm, sol_nw)
#            end
#        end
#    else
#        data["node_count"] = length(wm.data["junctions"]) + length(wm.data["reservoirs"])
#        data["link_count"] = length(wm.data["pipes"])
#
#        if status != :LocalInfeasible
#            solution_builder(wm, sol)
#        end
#    end
#
#    solution = Dict{String, Any}(
#        "solver" => string(typeof(wm.model.solver)),
#        "status" => status,
#        "objective" => objective,
#        "objective_lb" => guard_getobjbound(wm.model),
#        "solve_time" => solve_time,
#        "solution" => sol,
#        "machine" => Dict("cpu" => Sys.cpu_info()[1].model,
#                          "memory" => string(Sys.total_memory() / 2^30, " Gb")),
#        "data" => data
#    )
#
#    wm.solution = solution
#    return solution
#end
#
#function init_solution(wm::GenericWaterModel)
#    data_keys = ["multinetwork"]
#    return Dict{String,Any}(key => wm.data[key] for key in data_keys)
#end

function get_solution(wm::GenericWaterModel, sol::Dict{String,Any})
    add_setpoint(sol, wm, "pipes", "q", :q) # Get flow solutions.
    add_setpoint(sol, wm, "junctions", "h", :h) # Get head solutions (junctions).
    #add_setpoint(sol, wm, "reservoirs", "h", :h) # Get head solutions (reservoirs).
    return sol
end

function add_setpoint(sol, wm::GenericWaterModel, dict_name, param_name,
                      variable_symbol; index_name = "id",
                      default_value = (item) -> NaN, scale = (x, item) -> x,
                      extract_var = (var, idx, item) -> var[idx])
    sol_dict = get(sol, dict_name, Dict{String, Any}())

    if InfrastructureModels.ismultinetwork(wm.data)
        data_dict = wm.data["nw"]["$(wm.cnw)"][dict_name]
    else
        data_dict = wm.data[dict_name]
    end

    if length(data_dict) > 0
        sol[dict_name] = sol_dict
    end

    for (i, item) in data_dict
        idx = parse(Int, item[index_name])
        sol_item = sol_dict[i] = get(sol_dict, i, Dict{String, Any}())
        sol_item[param_name] = default_value(item)

        try
            # TODO: Use extract_var function to get the variable...
            variable = wm.var[:nw][wm.cnw][variable_symbol][idx]
            sol_item[param_name] = scale(getvalue(variable), item)
        catch
            # TODO: Put something here in case try fails...
        end
    end
end

#""
#function guard_objective_bound(model)
#    return JuMP.objective_bound(model)
#    try
#        return JuMP.objective_bound(model)
#    catch
#        if model.sense == MOI.SENSE_MIN
#            return -Inf
#        elseif model.sense == MOI.SENSE_MAX
#            return Inf
#        end
#    end
#end
