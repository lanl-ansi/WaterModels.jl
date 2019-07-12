""
function build_solution(wm::GenericWaterModel, solve_time; objective = NaN, solution_builder = get_solution)
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
    memory = string(Sys.total_memory() / 2^30, " GB")

    solution = Dict{String,Any}(
        "termination_status" => JuMP.termination_status(wm.model),
        "primal_status" => JuMP.primal_status(wm.model),
        "objective_value" => JuMP.objective_value(wm.model),
        "solve_time" => solve_time,
        "solution" => sol,
        "machine" => Dict("cpu" => cpu, "memory" => memory),
        "data" => data)
    #"dual_status" => JuMP.dual_status(wm.model),

    wm.solution = solution

    return solution
end

""
function init_solution(wm::GenericWaterModel)
    data_keys = ["per_unit"]
    return Dict{String,Any}(key => wm.data[key] for key in data_keys)
end

""
function get_solution(wm::GenericWaterModel, sol::Dict{String,<:Any})
    add_pipe_flow_rate_setpoint(sol, wm)
    add_pipe_resistance_setpoint(sol, wm)
    add_junction_head_setpoint(sol, wm)
end

function add_pipe_flow_rate_setpoint(sol, wm::GenericWaterModel)
    add_setpoint(sol, wm, "pipes", "q", :q)
end

function add_junction_head_setpoint(sol, wm::GenericWaterModel)
    if :h in keys(wm.var[:nw][wm.cnw])
        add_setpoint(sol, wm, "junctions", "h", :h)
    else
        # TODO: Add something here to compute heads from a flow rate solution.
    end
end

function add_pipe_resistance_setpoint(sol, wm::GenericWaterModel)
    if InfrastructureModels.ismultinetwork(wm.data)
        data_dict = wm.data["nw"]["$(wm.cnw)"]["pipes"]
    else
        data_dict = wm.data["pipes"]
    end

    sol_dict = get(sol, "pipes", Dict{String, Any}())

    if length(data_dict) > 0
        sol["pipes"] = sol_dict
    end

    if :x_res in keys(wm.var[:nw][wm.cnw])
        for (i, link) in data_dict
            a = link["id"]
            sol_item = sol_dict[i] = get(sol_dict, i, Dict{String, Any}())

            if a in keys(wm.var[:nw][wm.cnw][:x_res])
                x_res, r_id = findmax(JuMP.value.(wm.var[:nw][wm.cnw][:x_res][a]))
                sol_item["r"] = wm.ref[:nw][wm.cnw][:resistance][a][r_id]
            else
                x_res, r_id = findmin(wm.ref[:nw][wm.cnw][:resistance][a])
                sol_item["r"] = wm.ref[:nw][wm.cnw][:resistance][a][r_id]
            end
        end
    else
        for (i, link) in data_dict
            a = link["id"]
            sol_item = sol_dict[i] = get(sol_dict, i, Dict{String, Any}())
            x_res, r_id = findmin(wm.ref[:nw][wm.cnw][:resistance][a])
            sol_item["r"] = wm.ref[:nw][wm.cnw][:resistance][a][r_id]
        end
    end
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
        idx = item[index_name]
        sol_item = sol_dict[i] = get(sol_dict, i, Dict{String, Any}())
        sol_item[param_name] = default_value(item)

        try
            variable = extract_var(var(wm, wm.cnw, variable_symbol), idx, item)
            sol_item[param_name] = scale(JuMP.value(variable), item)
        catch
        end

        sol_item[param_name] = sol_item[param_name][1]
    end
end
