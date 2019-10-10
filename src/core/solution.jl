""
function build_solution(wm::AbstractWaterModel, solve_time; solution_builder=solution_owf!)
    sol = init_solution(wm)
    data = Dict{String, Any}("name" => wm.data["name"])
    valid_statuses = [_MOI.FEASIBLE_POINT, _MOI.NEARLY_FEASIBLE_POINT]

    if InfrastructureModels.ismultinetwork(wm.data)
        sol["multinetwork"] = true
        sol_nws = sol["nw"] = Dict{String,Any}()
        data_nws = data["nw"] = Dict{String,Any}()

        for (n, nw_data) in wm.data["nw"]
            sol_nw = sol_nws[n] = Dict{String,Any}()
            wm.cnw = parse(Int, n)

            if JuMP.primal_status(wm.model) in valid_statuses
                solution_builder(wm, sol_nw)
            end

            data_nws[n] = Dict("name" => get(nw_data, "name", "anonymous"),
                               "link_count" => length(ref(wm, :link)),
                               "node_count" => length(ref(wm, :node)))
        end
    else
        if JuMP.primal_status(wm.model) in valid_statuses
            solution_builder(wm, sol)
        end

        data["link_count"] = length(ref(wm, :link))
        data["node_count"] = length(ref(wm, :node))
    end

    solution = Dict{String,Any}(
        "optimizer" => JuMP.solver_name(wm.model),
        "termination_status" => JuMP.termination_status(wm.model),
        "primal_status" => JuMP.primal_status(wm.model),
        "dual_status" => JuMP.dual_status(wm.model),
        "objective" => _guard_objective_value(wm.model),
        "objective_lb" => _guard_objective_bound(wm.model),
        "solve_time" => solve_time,
        "solution" => sol,
        "machine" => Dict("cpu" => Sys.cpu_info()[1].model,
            "memory" => string(Sys.total_memory()/2^30, " Gb")),
        "data" => data)

    return solution
end

""
function init_solution(wm::AbstractWaterModel)
    data_keys = ["per_unit"]
    return Dict{String,Any}(key => wm.data[key] for key in data_keys)
end

function solution_owf!(wm::AbstractWaterModel, sol::Dict{String,<:Any})
    add_setpoint_node_head!(sol, wm)
    add_setpoint_pipe_flow!(sol, wm)
    add_setpoint_pipe_resistance!(sol, wm)
    add_setpoint_pump_flow!(sol, wm)
    add_setpoint_pump_status!(sol, wm)
    add_setpoint_reservoir!(sol, wm)
    add_setpoint_tank!(sol, wm)
end

""
function add_setpoint_node_head!(sol, wm::AbstractWaterModel)
    add_setpoint!(sol, wm, "node", "h", :h)
end

""
function add_setpoint_node_head!(sol, wm::AbstractCNLPModel)
    add_dual!(sol, wm, "node", "h", :flow_conservation, scale=(x, item) -> -x)
end

""
function add_setpoint_pipe_flow!(sol, wm::AbstractWaterModel)
    add_setpoint!(sol, wm, "pipe", "q", :q)
end

""
function add_setpoint_pipe_resistance!(sol, wm::AbstractWaterModel)
    if InfrastructureModels.ismultinetwork(wm.data)
        data_dict = wm.data["nw"]["$(wm.cnw)"]["pipe"]
    else
        data_dict = wm.data["pipe"]
    end

    sol_dict = get(sol, "pipe", Dict{String, Any}())

    if length(data_dict) > 0
        sol["pipe"] = sol_dict
    end

    if :x_res in keys(var(wm))
        for (i, link) in data_dict
            a = link["index"]
            sol_item = sol_dict[i] = get(sol_dict, i, Dict{String, Any}())

            if a in keys(var(wm, :x_res))
                x_res, r_id = findmax(JuMP.value.(var(wm, :x_res, a)))
                sol_item["r"] = ref(wm, :resistance, a)[r_id]
            else
                x_res, r_id = findmin(ref(wm, :resistance, a))
                sol_item["r"] = ref(wm, :resistance, a)[r_id]
            end
        end
    else
        for (i, link) in data_dict
            a = link["index"]
            sol_item = sol_dict[i] = get(sol_dict, i, Dict{String, Any}())
            x_res, r_id = findmin(ref(wm, :resistance, a))
            sol_item["r"] = ref(wm, :resistance, a)[r_id]
        end
    end
end

""
function add_setpoint_pump_flow!(sol, wm::AbstractWaterModel)
    add_setpoint!(sol, wm, "pump", "q", :q)
end

function add_setpoint_pump_status!(sol, wm::AbstractWaterModel)
    add_setpoint!(sol, wm, "pump", "x_pump", :x_pump)
end

function add_setpoint_reservoir!(sol, wm::AbstractWaterModel)
    add_setpoint!(sol, wm, "reservoir", "qr", :qr)
end

function add_setpoint_tank!(sol, wm::AbstractWaterModel)
    add_setpoint!(sol, wm, "tank", "qt", :qt)
end

"Adds values based on JuMP variables."
function add_setpoint!(
    sol,
    wm::AbstractWaterModel,
    dict_name,
    param_name,
    variable_symbol;
    index_name = "index",
    default_value = (item) -> NaN,
    scale = (x,item) -> x,
    var_key = (idx,item) -> idx,
    sol_dict = get(sol, dict_name, Dict{String,Any}()),
    status_name = "status",
    inactive_status_value = 0)

    has_variable_symbol = haskey(var(wm, wm.cnw), variable_symbol)

    variables = []
    if has_variable_symbol
        variables = var(wm, wm.cnw, variable_symbol)
    end

    if !has_variable_symbol || (!isa(variables, JuMP.VariableRef) && length(variables) == 0)
        add_setpoint_fixed!(sol, wm, dict_name, param_name; index_name=index_name, default_value=default_value)
        return
    end

    if InfrastructureModels.ismultinetwork(wm.data)
        data_dict = wm.data["nw"]["$(wm.cnw)"][dict_name]
    else
        data_dict = wm.data[dict_name]
    end

    if length(data_dict) > 0
        sol[dict_name] = sol_dict
    end

    for (i, item) in data_dict
        idx = Int(item[index_name])
        sol_item = sol_dict[i] = get(sol_dict, i, Dict{String,Any}())
        sol_item[param_name] = default_value(item)

        if item[status_name] != inactive_status_value
            var_id = var_key(idx, item)
            variables = var(wm, wm.cnw, variable_symbol)
            sol_item[param_name] = scale(JuMP.value(variables[var_id]), item)
        end
    end
end

"""
Adds setpoint values based on a given default_value function. This
significantly improves performance in models where values are not defined
"""
function add_setpoint_fixed!(
    sol,
    wm::AbstractWaterModel,
    dict_name,
    param_name;
    index_name = "index",
    default_value = (item) -> NaN,
    sol_dict = get(sol, dict_name, Dict{String,Any}()))

    if InfrastructureModels.ismultinetwork(wm.data)
        data_dict = wm.data["nw"]["$(wm.cnw)"][dict_name]
    else
        data_dict = wm.data[dict_name]
    end

    if length(data_dict) > 0
        sol[dict_name] = sol_dict
    end

    for (i,item) in data_dict
        idx = Int(item[index_name])
        sol_item = sol_dict[i] = get(sol_dict, i, Dict{String,Any}())
        sol_item[param_name] = default_value(item)
    end
end

"""
    function add_dual!(
        sol::AbstractDict,
        wm::AbstractWaterModel,
        dict_name::AbstractString,
        param_name::AbstractString,
        con_symbol::Symbol;
        index_name::AbstractString = "index",
        default_value::Function = (item) -> NaN,
        scale::Function = (x,item) -> x,
        con_key::Function = (idx,item) -> idx,
    )
This function takes care of adding the values of dual variables to the solution Dict.
# Arguments
- `sol::AbstractDict`: The dict where the desired final details of the solution are stored;
- `wm::AbstractWaterModel`: The WaterModel which has been considered;
- `dict_name::AbstractString`: The particular class of items for the solution (e.g. branch, bus);
- `param_name::AbstractString`: The name associated to the dual variable;
- `con_symbol::Symbol`: the Symbol attached to the class of constraints;
- `index_name::AbstractString = "index"`: ;
- `default_value::Function = (item) -> NaN`: a function that assign to each item a default value, for missing data;
- `scale::Function = (x,item) -> x`: a function to rescale the values of the dual variables, if needed;
- `con_key::Function = (idx,item) -> idx`: a method to extract the actual dual variables.
- `status_name::AbstractString: the status field of the given component type`
- `inactive_status_value::Any: the value of the status field indicating an inactive component`
"""
function add_dual!(
    sol::AbstractDict,
    wm::AbstractWaterModel,
    dict_name::AbstractString,
    param_name::AbstractString,
    con_symbol::Symbol;
    index_name::AbstractString = "index",
    default_value::Function = (item) -> NaN,
    scale::Function = (x, item) -> x,
    con_key::Function = (idx, item) -> idx,
    status_name = "status",
    inactive_status_value = 0)
    sol_dict = get(sol, dict_name, Dict{String,Any}())
    constraints = []
    has_con_symbol = haskey(con(wm, wm.cnw), con_symbol)

    if has_con_symbol
        constraints = con(wm, wm.cnw, con_symbol)
    end

    if !has_con_symbol || (!isa(constraints, JuMP.ConstraintRef) && length(constraints) == 0)
        add_dual_fixed!(sol, wm, dict_name, param_name; index_name=index_name, default_value=default_value)
        return
    end

    if ismultinetwork(wm)
        data_dict = wm.data["nw"]["$(wm.cnw)"][dict_name]
    else
        data_dict = wm.data[dict_name]
    end

    if length(data_dict) > 0
        sol[dict_name] = sol_dict
    end

    for (i, item) in data_dict
        idx = Int(item[index_name])
        sol_item = sol_dict[i] = get(sol_dict, i, Dict{String,Any}())
        sol_item[param_name] = default_value(item)

        if item[status_name] != inactive_status_value
            con_id = con_key(idx, item)
            constraints = con(wm, wm.cnw, con_symbol)
            sol_item[param_name] = scale(JuMP.dual(constraints[con_id]), item)
        end
    end
end

function add_dual_fixed!(
    sol::AbstractDict,
    wm::AbstractWaterModel,
    dict_name::AbstractString,
    param_name::AbstractString;
    index_name::AbstractString = "index",
    default_value::Function = (item) -> NaN)
    sol_dict = get(sol, dict_name, Dict{String,Any}())
    if ismultinetwork(wm)
        data_dict = wm.data["nw"]["$(wm.cnw)"][dict_name]
    else
        data_dict = wm.data[dict_name]
    end

    if length(data_dict) > 0
        sol[dict_name] = sol_dict
    end

    for (i,item) in data_dict
        idx = Int(item[index_name])
        sol_item = sol_dict[i] = get(sol_dict, i, Dict{String,Any}())
        sol_item[param_name] = default_value(item)
        sol_item[param_name] = sol_item[param_name][1]
    end
end

"Returns an objective value."
function _guard_objective_value(model::JuMP.AbstractModel)
    try
        return JuMP.objective_value(model)
    catch
        return JuMP.objective_sense(model) == _MOI.MAX_SENSE ? -Inf : Inf
    end
end


""
function _guard_objective_bound(model::JuMP.AbstractModel)
    try
        return JuMP.objective_bound(model)
    catch
        return JuMP.objective_sense(model) == _MOI.MAX_SENSE ? Inf : -Inf
    end
end
