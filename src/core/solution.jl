
"Build a water solution"
function build_solution{T}(wm::GenericWaterModel{T}, status, solve_time; objective = NaN, solution_builder = get_solution)
    if status != :Error
        objective = getobjectivevalue(wm.model)
        status = solver_status_dict(Symbol(typeof(wm.model.solver).name.module), status)
    end

    sol = solution_builder(wm)

    solution = Dict{AbstractString, Any}(
        "solver" => string(typeof(wm.model.solver)),
        "status" => status,
        "objective" => objective,
        "objective_lb" => guard_getobjbound(wm.model),
        "solve_time" => solve_time,
        "solution" => sol,
        "machine" => Dict("cpu" => Sys.cpu_info()[1].model,
                          "memory" => string(Sys.total_memory()/2^30, " Gb")),
        #"data" => Dict("junction_count" => length(wm.ref[:junction]),
        #               "connection_count" => length(wm.ref[:connection]))
    )

    wm.solution = solution

    return solution
end

""
function init_solution(wm::GenericWaterModel)
    return Dict{String,Any}()
end

" Get all the solution values "
function get_solution{T}(wm::GenericWaterModel{T})
    sol = init_solution(wm)
    #add_junction_pressure_setpoint(sol, wm)
    #add_connection_flow_setpoint(sol, wm)
    return sol
end

" Get the pressure squared solutions "
function add_junction_pressure_setpoint{T}(sol, wm::GenericWaterModel{T})
    add_setpoint(sol, wm, "junction", "p", :p; scale = (x,item) -> sqrt(x))
end

" Get the pressure squared solutions "
function add_load_setpoint{T}(sol, wm::GenericWaterModel{T})
    add_setpoint(sol, wm, "junction", "ql", :ql; default_value = (item) -> 0)
end

" Get the production set point "
function add_production_setpoint{T}(sol, wm::GenericWaterModel{T})
    add_setpoint(sol, wm, "junction", "qg", :qg; default_value = (item) -> 0)
end

" Get the direction set points"
function add_direction_setpoint{T}(sol, wm::GenericWaterModel{T})
    add_setpoint(sol, wm, "connection", "yp", :yp)
    add_setpoint(sol, wm, "connection", "yn", :yn)
end

" Get the valve solutions "
function add_valve_setpoint{T}(sol, wm::GenericWaterModel{T})
    add_setpoint(sol, wm, "connection", "valve", :v)
end

" Add the flow solutions "
function add_connection_flow_setpoint{T}(sol, wm::GenericWaterModel{T})
    add_setpoint(sol, wm, "connection", "f", :f)
end


function add_setpoint{T}(sol, wm::GenericWaterModel{T}, dict_name, param_name, variable_symbol; index_name = nothing, default_value = (item) -> NaN, scale = (x,item) -> x, extract_var = (var,idx,item) -> var[idx])
    sol_dict = get(sol, dict_name, Dict{String,Any}())
    if length(wm.data[dict_name]) > 0
        sol[dict_name] = sol_dict
    end

    for (i,item) in wm.data[dict_name]
        idx = parse(Int64,i)
        if index_name != nothing
            idx = Int(item[index_name])
        end
        sol_item = sol_dict[i] = get(sol_dict, i, Dict{String,Any}())
        sol_item[param_name] = default_value(item)
        try
            var = extract_var(wm.var[variable_symbol], idx, item)
            sol_item[param_name] = scale(getvalue(var), item)
        catch
        end
    end
end

solver_status_lookup = Dict{Any, Dict{Symbol, Symbol}}(
    :Ipopt => Dict(:Optimal => :LocalOptimal, :Infeasible => :LocalInfeasible),
    :ConicNonlinearBridge => Dict(:Optimal => :LocalOptimal, :Infeasible => :LocalInfeasible),
    # note that AmplNLWriter.AmplNLSolver is the solver type of bonmin
    :AmplNLWriter => Dict(:Optimal => :LocalOptimal, :Infeasible => :LocalInfeasible)
)

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
