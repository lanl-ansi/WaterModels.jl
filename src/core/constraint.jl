########################################################################
# This file defines commonly-used constraints for water systems models.
########################################################################

function constraint_flow_conservation{T}(wm::GenericWaterModel{T}, i, n::Int = wm.cnw)
    # Add the flow conservation constraints for junction nodes.
    out_arcs = collect(keys(filter((id, pipe) -> pipe["node1"] == i, wm.ref[:nw][n][:pipes])))
    out_vars = Array{JuMP.Variable}([wm.var[:nw][n][:q][a] for a in out_arcs])
    in_arcs = collect(keys(filter((id, pipe) -> pipe["node2"] == i, wm.ref[:nw][n][:pipes])))
    in_vars = Array{JuMP.Variable}([wm.var[:nw][n][:q][a] for a in in_arcs])

    if haskey(wm.ref[:nw][n][:junctions], i)
        demand = wm.ref[:nw][n][:junctions][i]["demand"]
        @constraint(wm.model, sum(in_vars) - sum(out_vars) == demand)
    elseif haskey(wm.ref[:nw][n][:reservoirs], i)
        junctions = values(wm.ref[:nw][n][:junctions])
        sum_demand = sum(junction["demand"] for junction in junctions)
        @constraint(wm.model, sum(out_vars) - sum(in_vars) >= 0.0)
        @constraint(wm.model, sum(out_vars) - sum(in_vars) <= sum_demand)
    end
end

function constraint_no_good{T}(wm::GenericWaterModel{T}, n::Int = wm.cnw)
    yp_solution = getvalue(wm.var[:nw][n][:yp])
    one_vars = [wm.var[:nw][n][:yp][idx[1]] for idx in keys(yp_solution) if yp_solution[idx[1]] >= 1]
    zero_vars = [wm.var[:nw][n][:yp][idx[1]] for idx in keys(yp_solution) if yp_solution[idx[1]] <= 0]
    @constraint(wm.model, sum(zero_vars) - sum(one_vars) >= 1 - length(one_vars))
end

function constraint_hw_diameter_selection{T}(wm::GenericWaterModel{T}, a, n::Int = wm.cnw)
    # Get variables associated with the discrete diameter choices.
    diameter_vars = wm.var[:nw][n][:diameter][a]

    # Get the corresponding diameter measurements (meters).
    diameters = [key[1] for key in keys(diameter_vars)]

    # Add a constraint that says at most one diameter may be selected.
    @constraint(wm.model, sum(diameter_vars) == 1)

    # Create a linear expression comprising the discrete selections for lambda.
    pipe = wm.ref[:nw][n][:ne_pipe][a]
    lambdas = [calc_friction_factor_hw_ne(pipe, d) for d in diameters]
    aff = AffExpr(diameter_vars[:], lambdas, 0.0)

    # Create an auxiliary variable that corresponds to the selected lambda.
    aux = @variable(wm.model, lowerbound = minimum(lambdas), upperbound = maximum(lambdas), start = minimum(lambdas))
    @constraint(wm.model, aux == aff)
end

function constraint_degree_two{T}(wm::GenericWaterModel{T}, idx, n::Int = wm.cnw)
    first = nothing
    last = nothing

    for i in wm.ref[:nw][n][:junction_connections][idx]
        connection = wm.ref[:nw][n][:connection][i]

        if connection["node1"] == idx
            other = connection["node2"]
        else
            other = connection["node1"]
        end

        if first == nothing
            first = other
        elseif first != other
            if last != nothing && last != other
                error(string("Error: adding a degree 2 constraint to a node with degree > 2: Junction ", idx))
            end

            last = other
        end
    end

    yp_first = filter(i -> wm.ref[:nw][n][:connection][i]["node1"] == first, wm.ref[:nw][n][:junction_connections][idx])
    yn_first = filter(i -> wm.ref[:nw][n][:connection][i]["node2"] == first, wm.ref[:nw][n][:junction_connections][idx])
    yp_last  = filter(i -> wm.ref[:nw][n][:connection][i]["node2"] == last,  wm.ref[:nw][n][:junction_connections][idx])
    yn_last  = filter(i -> wm.ref[:nw][n][:connection][i]["node1"] == last,  wm.ref[:nw][n][:junction_connections][idx])

    yp = wm.var[:nw][n][:yp]
    yn = wm.var[:nw][n][:yn]

    i = idx

    if !haskey(wm.con[:nw][n], :conserve_flow1)
        wm.con[:nw][n][:conserve_flow1] = Dict{String, ConstraintRef}()
        wm.con[:nw][n][:conserve_flow2] = Dict{String, ConstraintRef}()
        wm.con[:nw][n][:conserve_flow3] = Dict{String, ConstraintRef}()
        wm.con[:nw][n][:conserve_flow4] = Dict{String, ConstraintRef}()
    end

    if length(yn_first) > 0 && length(yp_last) > 0
        for i1 in yn_first
            for i2 in yp_last
                wm.con[:nw][n][:conserve_flow1][i] = @constraint(wm.model, yn[i1] == yp[i2])
                wm.con[:nw][n][:conserve_flow2][i] = @constraint(wm.model, yp[i1] == yn[i2])
                wm.con[:nw][n][:conserve_flow3][i] = @constraint(wm.model, yn[i1] + yn[i2] == 1)
                wm.con[:nw][n][:conserve_flow4][i] = @constraint(wm.model, yp[i1] + yp[i2] == 1)
            end
        end
    end

    if length(yn_first) > 0 && length(yn_last) > 0
        for i1 in yn_first
            for i2 in yn_last
                wm.con[:nw][n][:conserve_flow1][i] = @constraint(wm.model, yn[i1] == yn[i2])
                wm.con[:nw][n][:conserve_flow2][i] = @constraint(wm.model, yp[i1] == yp[i2])
                wm.con[:nw][n][:conserve_flow3][i] = @constraint(wm.model, yn[i1] + yp[i2] == 1)
                wm.con[:nw][n][:conserve_flow4][i] = @constraint(wm.model, yp[i1] + yn[i2] == 1)
            end
        end
    end

    if length(yp_first) > 0 && length(yp_last) > 0
        for i1 in yp_first
            for i2 in yp_last
                wm.con[:nw][n][:conserve_flow1][i] = @constraint(wm.model, yp[i1] == yp[i2])
                wm.con[:nw][n][:conserve_flow2][i] = @constraint(wm.model, yn[i1] == yn[i2])
                wm.con[:nw][n][:conserve_flow3][i] = @constraint(wm.model, yp[i1] + yn[i2] == 1)
                wm.con[:nw][n][:conserve_flow4][i] = @constraint(wm.model, yn[i1] + yp[i2] == 1)
            end
        end
    end

    if length(yp_first) > 0 && length(yn_last) > 0
        for i1 in yp_first
            for i2 in yn_last
                wm.con[:nw][n][:conserve_flow1][i] = @constraint(wm.model, yp[i1] == yn[i2])
                wm.con[:nw][n][:conserve_flow2][i] = @constraint(wm.model, yn[i1] == yp[i2])
                wm.con[:nw][n][:conserve_flow3][i] = @constraint(wm.model, yp[i1] + yp[i2] == 1)
                wm.con[:nw][n][:conserve_flow4][i] = @constraint(wm.model, yn[i1] + yn[i2] == 1)
            end
        end
    end
end
