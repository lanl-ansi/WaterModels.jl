# Define MILP (mixed-integer linear, relaxed program) implementations of water distribution models.

function variable_flow_piecewise_weights(wm::AbstractMILPModel, n::Int=wm.cnw)
    # Set the number of breakpoints used in each outer-approximation.
    num_breakpoints = :num_breakpoints in keys(wm.ext) ? wm.ext[:num_breakpoints] : 1

    var(wm, n)[:lambda] = JuMP.@variable(wm.model, [a in ids(wm, n, :link_fixed),
        k in 1:num_breakpoints], base_name="lambda[$(n)]", lower_bound=0.0,
        upper_bound=1.0, start=get_start(ref(wm, n, :link), a, "lambda_start", 0.0))
end

function variable_flow_piecewise_weights_ne(wm::AbstractMILPModel, n::Int=wm.cnw)
    # Set the number of breakpoints used in each outer-approximation.
    num_breakpoints = :num_breakpoints in keys(wm.ext) ? wm.ext[:num_breakpoints] : 1
    n_r = Dict(a => length(ref(wm, n, :resistance, a)) for a in ids(wm, n, :link_ne))

    var(wm, n)[:lambda] = JuMP.@variable(wm.model, [a in ids(wm, n, :link_ne),
        r in 1:n_r[a], k in 1:num_breakpoints],
        base_name="lambda[$(n)]", lower_bound=0.0, upper_bound=1.0,
        start=get_start(ref(wm, n, :link), a, "lambda_start", 0.0))
end

function variable_flow_piecewise_adjacency(wm::AbstractMILPModel, n::Int=wm.cnw)
    # Set the number of breakpoints used in each outer-approximation.
    num_breakpoints = :num_breakpoints in keys(wm.ext) ? wm.ext[:num_breakpoints] : 1
    var(wm, n)[:x_pw] = JuMP.@variable(wm.model, [a in ids(wm, n, :link_fixed),
        k in 1:num_breakpoints-1], base_name="x_pw[$(n)]", binary=true,
        start=get_start(ref(wm, n, :link_fixed), a, "x_pw_start", 0.0))
end

function variable_flow_piecewise_adjacency_ne(wm::AbstractMILPModel, n::Int=wm.cnw)
    # Set the number of breakpoints used in each outer-approximation.
    num_breakpoints = :num_breakpoints in keys(wm.ext) ? wm.ext[:num_breakpoints] : 1
    var(wm, n)[:x_pw] = JuMP.@variable(wm.model, [a in ids(wm, n, :link_ne),
        k in 1:num_breakpoints-1], base_name="x_pw[$(n)]", binary=true,
        start=get_start(ref(wm, n, :link_ne), a, "x_pw_start", 0.0))
end

function variable_flow(wm::AbstractMILPModel, n::Int=wm.cnw; bounded::Bool=true)
    # Create undirected flow variables (i.e., q).
    bounded ? variable_flow_bounded(wm, n) : variable_flow_unbounded(wm, n)

    # Create variables required for convex combination piecewise approximation.
    variable_flow_piecewise_weights(wm, n)
    variable_flow_piecewise_adjacency(wm, n)
end

"Create network expansion flow variables for undirected flow formulations."
function variable_flow_ne(wm::AbstractMILPModel, n::Int=wm.cnw; bounded::Bool=true)
    # Create undirected flow variables (i.e., q).
    bounded ? variable_flow_bounded_ne(wm, n) : variable_flow_unbounded_ne(wm, n)

    # Create expressions capturing the relationships among q, and q_ne.
    var(wm, n)[:q] = JuMP.@expression(wm.model, [a in ids(wm, n, :link_ne)],
        sum(var(wm, n, :q_ne, a)))

    # Create resistance binary variables.
    variable_resistance(wm, n)

    # Create variables required for convex combination piecewise approximation.
    variable_flow_piecewise_weights_ne(wm, n)
    variable_flow_piecewise_adjacency_ne(wm, n)
end

"Pump head gain constraint when the pump status is ambiguous."
function constraint_head_gain_pump(wm::AbstractMILPModel, n::Int, a::Int, node_fr::Int, node_to::Int, curve_fun::Array{Float64})
    # Fix reverse flow variable to zero (since this is a pump).
    q = var(wm, n, :q, a)

    g = var(wm, n, :g, a)
    h_i = var(wm, n, :h, node_fr)
    h_j = var(wm, n, :h, node_to)
    x_pump = var(wm, n, :x_pump, a)

    # Define the (relaxed) head gain caused by the pump.
    c_1 = JuMP.@NLconstraint(wm.model, curve_fun[1]*qp^2 + curve_fun[2]*qp + curve_fun[3]*x_pump <= g)

    # If the pump is off, decouple the head difference relationship.
    c_2 = JuMP.@constraint(wm.model, h_j - h_i - g <= 1.0e6 * (1 - x_pump))
    c_3 = JuMP.@constraint(wm.model, h_j - h_i - g >= -1.0e6 * (1 - x_pump))

    # If the pump is off, the flow along the pump must be zero.
    c_4 = JuMP.@constraint(wm.model, qp <= qp_ub * x_pump)
    c_5 = JuMP.@constraint(wm.model, qp >= 1.0e-6 * x_pump)

    # Append the constraint array.
    con(wm, n, :head_gain)[a] = [c_1, c_2, c_3, c_4, c_5]
end

function get_breakpoint_values(breakpoints::Array{Float64}, alpha::Float64)
    return [sign(x) * abs(x)^alpha for x in breakpoints]
end

function constraint_head_loss_pipe_ne(wm::AbstractMILPModel, n::Int, a::Int, alpha::Float64, node_fr::Int, node_to::Int, L::Float64, resistances)
    # Set the number of breakpoints used in each outer-approximation.
    n_b = :num_breakpoints in keys(wm.ext) ? wm.ext[:num_breakpoints] : 1

    if n_b > 0
        h_i = var(wm, n, :h, node_fr)
        h_j = var(wm, n, :h, node_to)
        lambda = var(wm, n, :lambda)
        x_pw = var(wm, n, :x_pw)
        lhs = JuMP.AffExpr(0.0)

        for (r_id, r) in enumerate(resistances)
            x_res = var(wm, n, :x_res, a)[r_id]
            q_ne = var(wm, n, :q_ne, a)[r_id]

            c_1 = JuMP.@constraint(wm.model, sum(lambda[a, r_id, k] for k in 1:n_b) == x_res)
            c_2 = JuMP.@constraint(wm.model, lambda[a, r_id, 1] <= x_pw[a, 1])
            c_3 = JuMP.@constraint(wm.model, lambda[a, r_id, n_b] <= x_pw[a, n_b-1])
            append!(con(wm, n, :head_loss, a), [c_1, c_2, c_3])

            for k in 2:n_b-1
                c_4 = JuMP.@constraint(wm.model, lambda[a, r_id, k] <= x_pw[a, k-1] + x_pw[a, k])
                append!(con(wm, n, :head_loss, a), [c_4])
            end

            q_ne = var(wm, n, :q_ne, a)[r_id]
            q_ne_lb = JuMP.has_lower_bound(q_ne) ? JuMP.lower_bound(q_ne) : -10.0
            q_ne_ub = JuMP.has_upper_bound(q_ne) ? JuMP.upper_bound(q_ne) : 10.0

            breakpoints = range(q_ne_lb, stop=q_ne_ub, length=n_b)
            f = get_breakpoint_values(collect(breakpoints), alpha)
            JuMP.add_to_expression!(lhs, r * sum(f[k] * lambda[a, r_id, k] for k in 1:n_b))

            q_ne_lhs = sum(breakpoints[k] * lambda[a, r_id, k] for k in 1:n_b)
            c_5 = JuMP.@constraint(wm.model, q_ne_lhs == q_ne)
            append!(con(wm, n, :head_loss, a), [c_5])

            # TODO: Address this hack for initializing piecewise binary variables.
            q_ne_start = JuMP.start_value(q_ne)

            if q_ne_start != 0.0
                for k in 1:n_b-1
                    if q_ne_start >= breakpoints[k] && q_ne_start <= breakpoints[k+1]
                        start_k = (q_ne_start - breakpoints[k]) / (breakpoints[k+1] - breakpoints[k])
                        start_kp1 = (breakpoints[k+1] - q_ne_start) / (breakpoints[k+1] - breakpoints[k])
                        JuMP.set_start_value(lambda[a, r_id, k], start_k)
                        JuMP.set_start_value(lambda[a, r_id, k+1], start_kp1)
                        JuMP.set_start_value(x_pw[a, k], 1.0)
                    end
                end

                if q_ne_start >= breakpoints[n_b]
                    JuMP.set_start_value(lambda[a, r_id, n_b], 1.0)
                    JuMP.set_start_value(x_pw[a, n_b-1], 1.0)
                end
            end
        end

        c_6 = JuMP.@constraint(wm.model, sum(x_pw[a, :]) == 1.0)
        c_7 = JuMP.@constraint(wm.model, lhs == inv(L) * (h_i - h_j))
        append!(con(wm, n, :head_loss, a), [c_6, c_7])
    end
end

function constraint_head_loss_pipe(wm::AbstractMILPModel, n::Int, a::Int, alpha::Float64, node_fr::Int, node_to::Int, L::Float64, r::Float64)
    # Set the number of breakpoints used in each outer-approximation.
    num_breakpoints = :num_breakpoints in keys(wm.ext) ? wm.ext[:num_breakpoints] : 1

    if num_breakpoints > 0
        q = var(wm, n, :q, a)
        h_i = var(wm, n, :h, node_fr)
        h_j = var(wm, n, :h, node_to)
        lambda = var(wm, n, :lambda)
        x_pw = var(wm, n, :x_pw)

        c_1 = JuMP.@constraint(wm.model, sum(lambda[a, :]) == 1.0)
        c_2 = JuMP.@constraint(wm.model, sum(x_pw[a, :]) == 1.0)
        c_3 = JuMP.@constraint(wm.model, lambda[a, 1] <= x_pw[a, 1])
        c_4 = JuMP.@constraint(wm.model, lambda[a, end] <= x_pw[a, end])

        q_lb = JuMP.has_lower_bound(q) ? JuMP.lower_bound(q) : -10.0
        q_ub = JuMP.has_upper_bound(q) ? JuMP.upper_bound(q) : 10.0
        breakpoints = range(q_lb, stop=q_ub, length=num_breakpoints)
        f = get_breakpoint_values(collect(breakpoints), alpha)

        lhs = r * sum(f[k] .* lambda[a, k] for k in 1:num_breakpoints)
        c_5 = JuMP.@constraint(wm.model, lhs == inv(L) * (h_i - h_j))

        # Append the constraint array.
        append!(con(wm, n, :head_loss, a), [c_1, c_2, c_3, c_4, c_5])

        for k in 2:num_breakpoints-1
            c_6 = JuMP.@constraint(wm.model, lambda[a, k] <= x_pw[a, k-1] + x_pw[a, k])
            append!(con(wm, n, :head_loss, a), [c_6])
        end
    end
end

function objective_wf(wm::AbstractMILPModel)
    JuMP.set_objective_sense(wm.model, _MOI.FEASIBILITY_SENSE)
end

function objective_owf(wm::AbstractMILPModel) 
    JuMP.set_objective_sense(wm.model, _MOI.FEASIBILITY_SENSE)
end
