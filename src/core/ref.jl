"""
Builds the ref dictionary from the data dictionary. Additionally the ref dictionary would
contain fields populated by the optional vector of ref_extensions provided as a keyword
argument.
"""
function build_ref(
    data::Dict{String,<:Any};
    ref_extensions::Vector{<:Function} = Vector{Function}([]),
)
    return _IM.build_ref(
        data,
        ref_add_core!,
        _wm_global_keys,
        ref_extensions = ref_extensions,
    )
end


function _calc_pump_energy_points(wm::AbstractWaterModel, nw::Int, pump_id::Int, num_points::Int)
    pump = ref(wm, nw, :pump, pump_id)
    constant = _DENSITY * _GRAVITY * ref(wm, nw, :time_step)
    curve_fun = _calc_head_curve_coefficients(pump)

    q_min, q_max = get(pump, "flow_min_forward", _FLOW_MIN), pump["flow_max"]
    q_build = range(q_min, stop = q_max, length = 100)
    f_build = _calc_cubic_flow_values(collect(q_build), curve_fun)

    if haskey(pump, "efficiency_curve")
        eff_curve = pump["efficiency_curve"]
        eff = _calc_efficiencies(collect(q_build), eff_curve)
    else
        eff = pump["efficiency"]
    end

    return q_build, constant .* inv.(eff) .* f_build
end


function _calc_pump_energy_ua(wm::AbstractWaterModel, nw::Int, pump_id::Int, q::Array{Float64, 1})
    q_true, f_true = _calc_pump_energy_points(wm, nw, pump_id, 100)
    f_interp = Interpolations.LinearInterpolation(q_true, f_true).(q)

    for i in 2:length(q)
        slope = (f_interp[i] - f_interp[i-1]) * inv(q[i] - q[i-1])
        true_ids = filter(x -> q_true[x] >= q[i-1] && q_true[x] <= q[i], 1:length(q_true))
        f_est_s = f_interp[i-1] .+ (slope .* (q_true[true_ids] .- q[i-1]))
        est_err = max(0.0, maximum(f_est_s .- f_true[true_ids]))
        f_interp[i-1:i] .-= est_err
    end

    return f_interp
end


function _calc_pump_energy(wm::AbstractWaterModel, nw::Int, pump_id::Int, q::Float64)
    pump = ref(wm, nw, :pump, pump_id)
    constant = _DENSITY * _GRAVITY * ref(wm, nw, :time_step)
    pc = _calc_head_curve_coefficients(pump)

    if haskey(pump, "efficiency_curve")
        eff_curve = pump["efficiency_curve"]
        eff = _calc_efficiencies([q], eff_curve)[1]
    else
        eff = pump["efficiency"]
    end

    return constant * inv(eff) * (pc[1]*q^3 + pc[2]*q^2 + pc[3]*q)
end


function _calc_head_loss_values(points::Array{Float64}, alpha::Float64)
    return [sign(x) * abs(x)^alpha for x in points]
end


function _calc_pump_gain_values(points::Array{Float64}, curve_fun::Array{Float64})
    return [curve_fun[1]*x*x + curve_fun[2]*x + curve_fun[3] for x in points]
end


function _calc_cubic_flow_values(points::Array{Float64}, curve_fun::Array{Float64})
    return [curve_fun[1]*x^3 + curve_fun[2]*x^2 + curve_fun[3]*x for x in points]
end


function _calc_efficiencies(points::Array{Float64}, curve::Array{Tuple{Float64, Float64}})
    q, eff = [[x[1] for x in curve], [x[2] for x in curve]]
    return Interpolations.LinearInterpolation(q, eff,
        extrapolation_bc=Interpolations.Flat()).(points)
end


function _get_function_from_head_curve(head_curve::Array{Tuple{Float64, Float64}})
    LsqFit.@. func(x, p) = p[1]*x*x + p[2]*x + p[3]

    if length(head_curve) > 1
        fit = LsqFit.curve_fit(func, first.(head_curve), last.(head_curve), [0.0, 0.0, 0.0])
        return LsqFit.coef(fit)
    elseif length(head_curve) == 1
        new_points = [(0.0, 1.33 * head_curve[1][2]), (2.0 * head_curve[1][1], 0.0)]
        head_curve = vcat(new_points, head_curve)
        fit = LsqFit.curve_fit(func, first.(head_curve), last.(head_curve), [0.0, 0.0, 0.0])
        return LsqFit.coef(fit)
    end
end


function _get_maximum_head_gain(pump::Dict{String,<:Any}, use_best_efficiency::Bool)
    if use_best_efficiency
        c = _calc_pump_best_efficiency_head_gain_curve(pump)
    else
        c = _calc_pump_head_gain_curve(pump)
    end

    q_at_max = -c[2] * inv(2.0 * c[1]) > 0.0 ? -c[2] * inv(2.0 * c[1]) : 0.0
    return c[1] * q_at_max^2 + c[2] * q_at_max + c[3]
end