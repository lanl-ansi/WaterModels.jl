########################################################################
# This file defines commonly-used constraints for water systems models.
########################################################################

function constraint_flow_conservation{T}(wm::GenericWaterModel{T}, i, n::Int = wm.cnw)
    # Add the flow conservation constraints for junction nodes.
    out_arcs = collect(keys(filter((id, pipe) -> pipe["node1"] == i, wm.ref[:nw][n][:pipes])))
    out_vars = Array{JuMP.Variable}([wm.var[:nw][n][:q][a] for a in out_arcs])
    in_arcs = collect(keys(filter((id, pipe) -> pipe["node2"] == i, wm.ref[:nw][n][:pipes])))
    in_vars = Array{JuMP.Variable}([wm.var[:nw][n][:q][a] for a in in_arcs])

    # Demands assume original units of liters per second.
    if haskey(wm.ref[:nw][n][:junctions], i)
        demand = wm.ref[:nw][n][:demand][i]
        @constraint(wm.model, sum(in_vars) - sum(out_vars) == demand)
    elseif haskey(wm.ref[:nw][n][:reservoirs], i)
        sum_demand = sum(values(wm.ref[:nw][n][:demand]))
        @constraint(wm.model, sum(out_vars) - sum(in_vars) >= 0.0)
        @constraint(wm.model, sum(out_vars) - sum(in_vars) <= sum_demand)
    end
end

function constraint_potential_flow_coupling{T}(wm::GenericWaterModel{T}, i, n::Int = wm.cnw)
    # Get the outgoing arcs of the node.
    arcs_from = collect(keys(filter((id, pipe) -> pipe["node1"] == i, wm.ref[:nw][n][:pipes])))
    headloss_type = wm.ref[:nw][n][:options]["headloss"]

    for a in arcs_from
        # Add the potential-flow coupling constraint.
        gamma = wm.var[:nw][n][:gamma][a]
        q = wm.var[:nw][n][:q][a]
        lambda = wm.ref[:nw][n][:lambda][a]

        if headloss_type == "h-w"
            # If Hazen-Williams formulation, use a piecewise outer-approximation.
            mids = linspace(getlowerbound(q), getupperbound(q), 3)
            mids_f = [lambda * (x^2)^0.926 for x in mids]
            mids_df = [lambda * 1.852*x / (x^2)^0.074 for x in mids]

				mids_df[isnan.(mids_df)] = 0.0
				int_1 = (mids_f[2] - mids_f[1]) / (mids_df[1] - mids_df[2])
				int_2 = (mids_f[3] - mids_f[2]) / (mids_df[2] - mids_df[3])

				mids = [x for x in mids]
				sort!(append!(mids, [int_1, int_2]))
				mids_f = [lambda * (x^2)^0.926 for x in mids]

            rhs = piecewiselinear(wm.model, q, mids, mids_f, method = :ZigZagInteger)
            @constraint(wm.model, gamma == 1.00 * rhs)
        elseif headloss_type == "d-w"
            # If Darcy-Weisbach formulation, we can probably handle it natively.
            @NLconstraint(wm.model, gamma >= lambda * q^2)
        end
    end
end

function constraint_define_gamma{T}(wm::GenericWaterModel{T}, a, n::Int = wm.cnw)
    # Get source and target nodes corresponding to the arc.
    i = wm.ref[:nw][n][:pipes][a]["node1"]
    j = wm.ref[:nw][n][:pipes][a]["node2"]

    # Collect variables needed for the constraint.
    gamma = wm.var[:nw][n][:gamma][a]
    y_p = wm.var[:nw][n][:yp][a]
    y_n = wm.var[:nw][n][:yn][a]
    h_i = wm.var[:nw][n][:h][i]
    h_j = wm.var[:nw][n][:h][j]

    # Get additional data related to the variables.
    h_i_lb = getlowerbound(h_i)
    h_i_ub = getupperbound(h_i)
    h_j_lb = getlowerbound(h_j)
    h_j_ub = getupperbound(h_j)

    # Add the required constraints.
    @constraint(wm.model, h_j - h_i + (h_i_lb - h_j_ub) * (y_p - y_n + 1) <= gamma)
    @constraint(wm.model, h_i - h_j + (h_i_ub - h_j_lb) * (y_p - y_n - 1) <= gamma)
    @constraint(wm.model, h_j - h_i + (h_i_ub - h_j_lb) * (y_p - y_n + 1) >= gamma)
    @constraint(wm.model, h_i - h_j + (h_i_lb - h_j_ub) * (y_p - y_n - 1) >= gamma)
end

function constraint_bidirectional_flow{T}(wm::GenericWaterModel{T}, a, n::Int = wm.cnw)
    # Get source and target nodes corresponding to the arc.
    i = wm.ref[:nw][n][:pipes][a]["node1"]
    j = wm.ref[:nw][n][:pipes][a]["node2"]

    # Collect variables needed for the constraint.
    q = wm.var[:nw][n][:q][a]
    y_p = wm.var[:nw][n][:yp][a]
    y_n = wm.var[:nw][n][:yn][a]
    h_i = wm.var[:nw][n][:h][i]
    h_j = wm.var[:nw][n][:h][j]

    # Get additional data related to the variables.
    h_i_lb = getlowerbound(h_i)
    h_i_ub = getupperbound(h_i)
    h_j_lb = getlowerbound(h_j)
    h_j_ub = getupperbound(h_j)

    # Add the first set of constraints.
    sum_demand = sum(values(wm.ref[:nw][n][:demand]))
    @constraint(wm.model, (y_p - 1) * sum_demand <= q)
    @constraint(wm.model, (1 - y_n) * sum_demand >= q)

    # Add the second set of constraints.
    @constraint(wm.model, (1 - y_p) * (h_i_lb - h_j_ub) <= h_i - h_j)
    @constraint(wm.model, (1 - y_n) * (h_i_ub - h_j_lb) >= h_i - h_j)

    # Add the third set of constraints.
    @constraint(wm.model, y_p + y_n == 1)
end
