# This file is added as a temporary fix to work with scip.
# In future, the pipe head_loss function needs to be changed

function constraint_pipe_head_loss_scip_lp(
    wm::AbstractWaterModel,
    a::Int;
    nw::Int = nw_id_default,
    kwargs...,
)
    node_fr = ref(wm, nw, :pipe, a, "node_fr")
    node_to = ref(wm, nw, :pipe, a, "node_to")
    exponent = _get_exponent_from_head_loss_form(wm.ref[:it][wm_it_sym][:head_loss])
    L = ref(wm, nw, :pipe, a, "length")

    wm_data = get_wm_data(wm.data)
    head_loss, viscosity = wm_data["head_loss"], wm_data["viscosity"]
    base_length = wm_data["per_unit"] ? wm_data["base_length"] : 1.0
    base_mass = wm_data["per_unit"] ? wm_data["base_mass"] : 1.0
    base_time = wm_data["per_unit"] ? wm_data["base_time"] : 1.0

    r = _calc_pipe_resistance(
        ref(wm, nw, :pipe, a),
        head_loss,
        viscosity,
        base_length,
        base_mass,
        base_time,
    )
    q_max_reverse = min(get(ref(wm, nw, :pipe, a), "flow_max_reverse", 0.0), 0.0)
    q_min_forward = max(get(ref(wm, nw, :pipe, a), "flow_min_forward", 0.0), 0.0)

    _initialize_con_dict(wm, :pipe_head_loss, nw = nw, is_array = true)
    con(wm, nw, :pipe_head_loss)[a] = Array{JuMP.ConstraintRef}([])
    constraint_pipe_head_loss_scip_lp(
        wm,
        nw,
        a,
        node_fr,
        node_to,
        exponent,
        L,
        r,
        q_max_reverse,
        q_min_forward,
    )
end



function constraint_pipe_head_loss_scip_lp(
    wm::AbstractNCDModel,
    n::Int,
    a::Int,
    node_fr::Int,
    node_to::Int,
    exponent::Float64,
    L::Float64,
    r::Float64,
    q_max_reverse::Float64,
    q_min_forward::Float64,
)
# Get the variable for flow directionality.
y = var(wm, n, :y_pipe, a)

# Get variables for positive flow and head difference.
qp, dhp = var(wm, n, :qp_pipe, a), var(wm, n, :dhp_pipe, a)

# Get positively-directed convex combination variables.
lambda_p, x_p = var(wm, n, :lambda_p_pipe), var(wm, n, :x_p_pipe)

# Get the corresponding positive flow partitioning.
partition_p = get_pipe_flow_partition_positive(ref(wm, n, :pipe, a))
bp_range, bp_range_m1 = 1:length(partition_p), 1:length(partition_p)-1

# Add constraints for the positive flow piecewise approximation.
c_1 = JuMP.@constraint(wm.model, sum(lambda_p[a, k] for k in bp_range) == y)
qp_sum = sum(partition_p[k] * lambda_p[a, k] for k in bp_range)
scalar = _get_scaling_factor(vcat(qp_sum.terms.vals, [1.0]))
c_2 = JuMP.@constraint(wm.model, scalar * qp_sum == scalar * qp)
append!(con(wm, n, :pipe_head_loss, a), [c_1, c_2])

if length(partition_p) > 1
    # If there are multiple points, constrain the convex combination.
    c_3 = JuMP.@constraint(wm.model, sum(x_p[a, k] for k in bp_range_m1) == y)
    c_4 = JuMP.@constraint(wm.model, lambda_p[a, 1] <= x_p[a, 1])
    c_5 = JuMP.@constraint(wm.model, lambda_p[a, bp_range[end]] <= x_p[a, bp_range_m1[end]])
    append!(con(wm, n, :pipe_head_loss, a), [c_3, c_4, c_5])
end

# Add a constraint that upper-bounds the head loss variable.
if maximum(partition_p) != 0.0
    f_p = r * partition_p.^exponent
    f_p_ub_expr = sum(f_p[k] * lambda_p[a, k] for k in bp_range)
    scalar = _get_scaling_factor(vcat(f_p_ub_expr.terms.vals, [1.0 / L]))
    c_6 = JuMP.@constraint(wm.model, scalar * dhp / L <= scalar * f_p_ub_expr)
    append!(con(wm, n, :pipe_head_loss, a), [c_6])
else
    c_6 = JuMP.@constraint(wm.model, dhp == 0.0)
    append!(con(wm, n, :pipe_head_loss, a), [c_6])
end

# Loop over consequential points (i.e., those that have nonzero head loss).
for flow_value in filter(x -> x > 0.0, partition_p)
    # Add head loss outer (i.e., lower) approximations.
    lhs_p = r * _calc_head_loss_oa(qp, y, flow_value, exponent)

    if minimum(abs.(lhs_p.terms.vals)) >= 1.0e-4
        scalar = _get_scaling_factor(vcat(lhs_p.terms.vals, [1.0 / L]))
        c_7 = JuMP.@constraint(wm.model, scalar * lhs_p <= scalar * dhp / L)
        append!(con(wm, n, :pipe_head_loss, a), [c_7])
    end
end

for k in 2:length(partition_p)-1
    # Add the adjacency constraints for piecewise variables.
    c_8 = JuMP.@constraint(wm.model, lambda_p[a, k] <= x_p[a, k-1] + x_p[a, k])
    append!(con(wm, n, :pipe_head_loss, a), [c_8])
end

# Get variables for negative flow and head difference.
qn, dhn = var(wm, n, :qn_pipe, a), var(wm, n, :dhn_pipe, a)

# Get negatively-directed convex combination variable.
lambda_n, x_n = var(wm, n, :lambda_n_pipe), var(wm, n, :x_n_pipe)

# Get the corresponding negative flow partitioning (negated).
partition_n = sort(-get_pipe_flow_partition_negative(ref(wm, n, :pipe, a)))
bn_range, bn_range_m1 = 1:length(partition_n), 1:length(partition_n)-1

# Add constraints for the negative flow piecewise approximation.
c_9 = JuMP.@constraint(wm.model, sum(lambda_n[a, k] for k in bn_range) == 1.0 - y)
qn_sum = sum(partition_n[k] * lambda_n[a, k] for k in bn_range)
scalar = _get_scaling_factor(vcat(qn_sum.terms.vals, [1.0]))
c_10 = JuMP.@constraint(wm.model, scalar * qn_sum == scalar * qn)
append!(con(wm, n, :pipe_head_loss, a), [c_9, c_10])

if length(partition_n) > 1
    # If there are multiple points, constrain the convex combination.
    c_11 = JuMP.@constraint(wm.model, sum(x_n[a, k] for k in bn_range_m1) == 1.0 - y)
    c_12 = JuMP.@constraint(wm.model, lambda_n[a, 1] <= x_n[a, 1])
    c_13 = JuMP.@constraint(wm.model, lambda_n[a, bn_range[end]] <= x_n[a, bn_range_m1[end]])
    append!(con(wm, n, :pipe_head_loss, a), [c_11, c_12, c_13])
end

# Add a constraint that upper-bounds the head loss variable.
if maximum(partition_n) != 0.0
    f_n = r .* partition_n.^exponent
    f_n_ub_expr = sum(f_n[k] * lambda_n[a, k] for k in bn_range)
    scalar = _get_scaling_factor(vcat(f_n_ub_expr.terms.vals, [1.0 / L]))
    c_14 = JuMP.@constraint(wm.model, scalar * dhn / L <= scalar * f_n_ub_expr)
    append!(con(wm, n, :pipe_head_loss, a), [c_14])
else
    c_14 = JuMP.@constraint(wm.model, dhn == 0.0)
    append!(con(wm, n, :pipe_head_loss, a), [c_14])
end

# Loop over consequential points (i.e., those that have nonzero head loss).
for flow_value in filter(x -> x > 0.0, partition_n)
    # Add head loss outer (i.e., lower) approximations.
    lhs_n = r * _calc_head_loss_oa(qn, 1.0 - y, flow_value, exponent)

    if minimum(abs.(lhs_n.terms.vals)) >= 1.0e-4
        scalar = _get_scaling_factor(vcat(lhs_n.terms.vals, [1.0 / L]))
        c_15 = JuMP.@constraint(wm.model, scalar * lhs_n <= scalar * dhn / L)
        append!(con(wm, n, :pipe_head_loss, a), [c_15])
    end
end

for k in 2:length(partition_n)-1
    # Add the adjacency constraints for piecewise variables.
    c_16 = JuMP.@constraint(wm.model, lambda_n[a, k] <= x_n[a, k-1] + x_n[a, k])
    append!(con(wm, n, :pipe_head_loss, a), [c_16])
end
end



# function variable_flow_scip(
#     wm::AbstractNCDModel;
#     nw::Int = nw_id_default,
#     bounded::Bool = true,
#     report::Bool = true,
# )
#     for name in _LINK_COMPONENTS
#         # Create directed flow (`qp` and `qn`) variables for each component.
#         _variable_component_flow(wm, name; nw = nw, bounded = bounded, report = report)
#
#         # Create directed flow binary direction variables (`y`) for each component.
#         _variable_component_direction(wm, name; nw = nw, report = report)
#     end
#
#
#         # Create variables involved in convex combination constraints for pipes.
#         pipe_partition_p = Dict{Int, Vector{Float64}}(a =>
#             get_pipe_flow_partition_positive(pipe)
#             for (a, pipe) in ref(wm, nw, :pipe))
#
#         var(wm, nw)[:lambda_p_pipe] = JuMP.@variable(wm.model,
#             [a in ids(wm, nw, :pipe), k in 1:length(pipe_partition_p[a])],
#             base_name = "$(nw)_lambda_p", lower_bound = 0.0, upper_bound = 1.0,
#             start = comp_start_value(ref(wm, nw, :pipe, a), "lambda_p_start", k))
#
#         var(wm, nw)[:x_p_pipe] = JuMP.@variable(wm.model,
#             [a in ids(wm, nw, :pipe), k in 1:length(pipe_partition_p[a])-1],
#             base_name = "$(nw)_x_p", binary = true,
#             start = comp_start_value(ref(wm, nw, :pipe, a), "x_p_start"))
#
#         pipe_partition_n = Dict{Int, Vector{Float64}}(a =>
#             sort(-get_pipe_flow_partition_negative(pipe))
#             for (a, pipe) in ref(wm, nw, :pipe))
#
#         var(wm, nw)[:lambda_n_pipe] = JuMP.@variable(wm.model,
#             [a in ids(wm, nw, :pipe), k in 1:length(pipe_partition_n[a])],
#             base_name = "$(nw)_lambda_n", lower_bound = 0.0, upper_bound = 1.0,
#             start = comp_start_value(ref(wm, nw, :pipe, a), "lambda_n_start", k))
#
#         var(wm, nw)[:x_n_pipe] = JuMP.@variable(wm.model,
#             [a in ids(wm, nw, :pipe), k in 1:length(pipe_partition_n[a])-1],
#             base_name = "$(nw)_x_n", binary = true,
#             start = comp_start_value(ref(wm, nw, :pipe, a), "x_n_start"))
#
#     for name in ["des_pipe", "pipe"]
#         # Create directed head difference (`dhp` and `dhn`) variables for each component.
#         _variable_component_head_difference(
#             wm,
#             name;
#             nw = nw,
#             bounded = bounded,
#             report = report,
#         )
#     end
# end

# function constraint_pipe_head_loss_scip(
#     wm::AbstractNCDModel,
#     n::Int,
#     a::Int,
#     node_fr::Int,
#     node_to::Int,
#     exponent::Float64,
#     L::Float64,
#     r::Float64,
#     q_max_reverse::Float64,
#     q_min_forward::Float64
# )
#     # Get the variable for flow directionality.
#     y = var(wm, n, :y_pipe, a)
#
#     # Get variables for positive flow and head difference.
#     qp, dhp = var(wm, n, :qp_pipe, a), var(wm, n, :dhp_pipe, a)
#
#     # Get positively-directed convex combination variables.
#     lambda_p, x_p = var(wm, n, :lambda_p_pipe), var(wm, n, :x_p_pipe)
#
#     # Get the corresponding positive flow partitioning.
#     partition_p = get_pipe_flow_partition_positive(ref(wm, n, :pipe, a))
#     bp_range, bp_range_m1 = 1:length(partition_p), 1:length(partition_p)-1
#
#     # Add constraints for the positive flow piecewise approximation.
#     c_1 = JuMP.@constraint(wm.model, sum(lambda_p[a, k] for k in bp_range) == y)
#     qp_sum = sum(partition_p[k] * lambda_p[a, k] for k in bp_range)
#     scalar = _get_scaling_factor(vcat(qp_sum.terms.vals, [1.0]))
#     c_2 = JuMP.@constraint(wm.model, scalar * qp_sum == scalar * qp)
#     append!(con(wm, n, :pipe_head_loss, a), [c_1, c_2])
#
#     if length(partition_p) > 1
#         # If there are multiple points, constrain the convex combination.
#         c_3 = JuMP.@constraint(wm.model, sum(x_p[a, k] for k in bp_range_m1) == y)
#         c_4 = JuMP.@constraint(wm.model, lambda_p[a, 1] <= x_p[a, 1])
#         c_5 = JuMP.@constraint(wm.model, lambda_p[a, bp_range[end]] <= x_p[a, bp_range_m1[end]])
#         append!(con(wm, n, :pipe_head_loss, a), [c_3, c_4, c_5])
#     end
#
#     # Add a constraint that upper-bounds the head loss variable.
#     if maximum(partition_p) != 0.0
#         f_p = r * partition_p.^exponent
#         f_p_ub_expr = sum(f_p[k] * lambda_p[a, k] for k in bp_range)
#         scalar = _get_scaling_factor(vcat(f_p_ub_expr.terms.vals, [1.0 / L]))
#         c_6 = JuMP.@constraint(wm.model, scalar * dhp / L <= scalar * f_p_ub_expr)
#         append!(con(wm, n, :pipe_head_loss, a), [c_6])
#     else
#         c_6 = JuMP.@constraint(wm.model, dhp == 0.0)
#         append!(con(wm, n, :pipe_head_loss, a), [c_6])
#     end
#
#     # Loop over consequential points (i.e., those that have nonzero head loss).
#     for flow_value in filter(x -> x > 0.0, partition_p)
#         # Add head loss outer (i.e., lower) approximations.
#         lhs_p = r * _calc_head_loss_oa(qp, y, flow_value, exponent)
#
#         if minimum(abs.(lhs_p.terms.vals)) >= 1.0e-4
#             scalar = _get_scaling_factor(vcat(lhs_p.terms.vals, [1.0 / L]))
#             c_7 = JuMP.@constraint(wm.model, scalar * lhs_p <= scalar * dhp / L)
#             append!(con(wm, n, :pipe_head_loss, a), [c_7])
#         end
#     end
#
#     for k in 2:length(partition_p)-1
#         # Add the adjacency constraints for piecewise variables.
#         c_8 = JuMP.@constraint(wm.model, lambda_p[a, k] <= x_p[a, k-1] + x_p[a, k])
#         append!(con(wm, n, :pipe_head_loss, a), [c_8])
#     end
#
#     # Get variables for negative flow and head difference.
#     qn, dhn = var(wm, n, :qn_pipe, a), var(wm, n, :dhn_pipe, a)
#
#     # Get negatively-directed convex combination variable.
#     lambda_n, x_n = var(wm, n, :lambda_n_pipe), var(wm, n, :x_n_pipe)
#
#     # Get the corresponding negative flow partitioning (negated).
#     partition_n = sort(-get_pipe_flow_partition_negative(ref(wm, n, :pipe, a)))
#     bn_range, bn_range_m1 = 1:length(partition_n), 1:length(partition_n)-1
#
#     # Add constraints for the negative flow piecewise approximation.
#     c_9 = JuMP.@constraint(wm.model, sum(lambda_n[a, k] for k in bn_range) == 1.0 - y)
#     qn_sum = sum(partition_n[k] * lambda_n[a, k] for k in bn_range)
#     scalar = _get_scaling_factor(vcat(qn_sum.terms.vals, [1.0]))
#     c_10 = JuMP.@constraint(wm.model, scalar * qn_sum == scalar * qn)
#     append!(con(wm, n, :pipe_head_loss, a), [c_9, c_10])
#
#     if length(partition_n) > 1
#         # If there are multiple points, constrain the convex combination.
#         c_11 = JuMP.@constraint(wm.model, sum(x_n[a, k] for k in bn_range_m1) == 1.0 - y)
#         c_12 = JuMP.@constraint(wm.model, lambda_n[a, 1] <= x_n[a, 1])
#         c_13 = JuMP.@constraint(wm.model, lambda_n[a, bn_range[end]] <= x_n[a, bn_range_m1[end]])
#         append!(con(wm, n, :pipe_head_loss, a), [c_11, c_12, c_13])
#     end
#
#     # Add a constraint that upper-bounds the head loss variable.
#     if maximum(partition_n) != 0.0
#         f_n = r .* partition_n.^exponent
#         f_n_ub_expr = sum(f_n[k] * lambda_n[a, k] for k in bn_range)
#         scalar = _get_scaling_factor(vcat(f_n_ub_expr.terms.vals, [1.0 / L]))
#         c_14 = JuMP.@constraint(wm.model, scalar * dhn / L <= scalar * f_n_ub_expr)
#         append!(con(wm, n, :pipe_head_loss, a), [c_14])
#     else
#         c_14 = JuMP.@constraint(wm.model, dhn == 0.0)
#         append!(con(wm, n, :pipe_head_loss, a), [c_14])
#     end
#
#     # Loop over consequential points (i.e., those that have nonzero head loss).
#     for flow_value in filter(x -> x > 0.0, partition_n)
#         # Add head loss outer (i.e., lower) approximations.
#         lhs_n = r * _calc_head_loss_oa(qn, 1.0 - y, flow_value, exponent)
#
#         if minimum(abs.(lhs_n.terms.vals)) >= 1.0e-4
#             scalar = _get_scaling_factor(vcat(lhs_n.terms.vals, [1.0 / L]))
#             c_15 = JuMP.@constraint(wm.model, scalar * lhs_n <= scalar * dhn / L)
#             append!(con(wm, n, :pipe_head_loss, a), [c_15])
#         end
#     end
#
#     for k in 2:length(partition_n)-1
#         # Add the adjacency constraints for piecewise variables.
#         c_16 = JuMP.@constraint(wm.model, lambda_n[a, k] <= x_n[a, k-1] + x_n[a, k])
#         append!(con(wm, n, :pipe_head_loss, a), [c_16])
#     end
# end
