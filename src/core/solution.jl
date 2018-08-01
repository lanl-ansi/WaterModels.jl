function build_solution{T}(wm::GenericWaterModel{T}, status, solve_time;
                           objective = NaN, solution_builder = get_solution)
    if status != :Error
        objective = getobjectivevalue(wm.model)
        status = solver_status_dict(Symbol(typeof(wm.model.solver).name.module), status)
    end

    sol = init_solution(wm)
    data = Dict{String, Any}("name" => wm.data["title"])

    if wm.data["multinetwork"]
        sol_nws = sol["nw"] = Dict{String, Any}()
        data_nws = data["nw"] = Dict{String, Any}()

        for (n, nw_data) in wm.data["nw"]
            sol_nw = sol_nws[n] = Dict{String, Any}()
            wm.cnw = parse(Int, n)
            data_nws[n] = Dict(
                "name" => get(nw_data, "name", "anonymous"),
                "node_count" => length(nw_data["junctions"]) + length(nw_data["reservoirs"]),
                "link_count" => length(nw_data["pipes"])
            )

            if status != :LocalInfeasible
                solution_builder(wm, sol_nw)
            end
        end
    else
        data["node_count"] = length(wm.data["junctions"]) + length(wm.data["reservoirs"])
        data["link_count"] = length(wm.data["pipes"])

        if status != :LocalInfeasible
            solution_builder(wm, sol)
        end
    end

    solution = Dict{String, Any}(
        "solver" => string(typeof(wm.model.solver)),
        "status" => status,
        "objective" => objective,
        "objective_lb" => guard_getobjbound(wm.model),
        "solve_time" => solve_time,
        "solution" => sol,
        "machine" => Dict("cpu" => Sys.cpu_info()[1].model,
                          "memory" => string(Sys.total_memory() / 2^30, " Gb")),
        "data" => data
    )

    wm.solution = solution
    return solution
end

function init_solution(wm::GenericWaterModel)
    data_keys = ["multinetwork"]
    return Dict{String,Any}(key => wm.data[key] for key in data_keys)
end

function get_solution(wm::GenericWaterModel, sol::Dict{String,Any})
    add_setpoint(sol, wm, "pipes", "gamma", :gamma) # Get absolute value of head difference.
    add_setpoint(sol, wm, "pipes", "q", :q) # Get flow.
    add_setpoint(sol, wm, "junctions", "h", :h) # Get head solution (junctions).
    add_setpoint(sol, wm, "reservoirs", "h", :h) # Get head solution (reservoirs).
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
        idx = String(item[index_name]) # TODO: idx needs to be an integer.
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

solver_status_lookup = Dict{Any, Dict{Symbol, Symbol}}(
    :AmplNLWriter => Dict(:Optimal => :LocalOptimal, :Infeasible => :LocalInfeasible),
    :ConicNonlinearBridge => Dict(:Optimal => :LocalOptimal, :Infeasible => :LocalInfeasible),
    :Gurobi => Dict(:Optimal => :LocalOptimal, :Infeasible => :LocalInfeasible),
    :Ipopt => Dict(:Optimal => :LocalOptimal, :Infeasible => :LocalInfeasible),
    :Pajarito => Dict(:Optimal => :LocalOptimal, :Infeasible => :LocalInfeasible))

# translates solver status codes to our status codes
function solver_status_dict(solver_type, status)
    for (st, solver_stat_dict) in solver_status_lookup
        if solver_type == st
            if status in keys(solver_stat_dict)
                return solver_stat_dict[status]
            else
                return status
            end
        end
    end

    return status
end

function guard_getobjbound(model)
    try
        getobjbound(model)
    catch
        -Inf
    end
end
