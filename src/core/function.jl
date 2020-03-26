###############################################################################
# This file defines the nonlinear head loss functions for water systems models.
###############################################################################

function if_alpha(alpha::Float64; convex::Bool=false)
    return function(x::Float64)
        return inv(2.0 + alpha) * (x*x)^(1.0 + 0.5*alpha)
    end
end

function f_alpha(alpha::Float64; convex::Bool=false)
    if convex
        return function(x::Float64)
            return (x*x)^(0.5 + 0.5*alpha)
        end
    else
        return function(x::Float64)
            return sign(x) * (x*x)^(0.5 + 0.5*alpha)
        end
    end
end

function df_alpha(alpha::Float64; convex::Bool=false)
    if convex
        return function(x::Float64)
            return (1.0 + alpha) * sign(x) * (x*x)^(0.5*alpha)
        end
    else
        return function(x::Float64)
            return (1.0 + alpha) * (x*x)^(0.5*alpha)
        end
    end
end

function d2f_alpha(alpha::Float64; convex::Bool=false)
    if convex
        return function(x::Float64)
            return alpha * (1.0 + alpha) * (x*x)^(0.5*alpha - 0.5)
        end
    else
        return function(x::Float64)
            return sign(x) * alpha * (1.0 + alpha) * (x*x)^(0.5*alpha - 0.5)
        end
    end
end

function f_dual(alpha::Float64)
    return function(x::Float64)
        return alpha * inv(1.0 + alpha) * (x*x)^(0.5 + 0.5*inv(alpha))
    end
end

function df_dual(alpha::Float64)
    return function(x::Float64)
        return sign(x) * (x*x)^(0.5*inv(alpha))
    end
end

function d2f_dual(alpha::Float64)
    return function(x::Float64)
        return x != 0.0 ? inv(alpha) * (x*x)^(0.5*inv(alpha) - 0.5) : prevfloat(Inf)
    end
end

function f_primal(alpha::Float64)
    return function(x::Float64)
        return inv(1.0 + alpha) * (x*x)^(0.5 + 0.5*alpha)
    end
end

function df_primal(alpha::Float64)
    return function(x::Float64)
        return sign(x) * (x*x)^(0.5*alpha)
    end
end

function d2f_primal(alpha::Float64)
    return function(x::Float64)
        return alpha * (x*x)^(0.5*alpha - 0.5)
    end
end

function get_alpha_min_1(wm::AbstractWaterModel)
    alpha = [ref(wm, nw, :alpha) for nw in nw_ids(wm)]

    if !all(y -> y == alpha[1], alpha)
        Memento.error(_LOGGER, "Head loss exponents are different across the multinetwork.")
    else
        return alpha[1] - 1.0
    end
end

function head_loss_args(wm::AbstractCNLPModel)
    alpha_m1 = get_alpha_min_1(wm)
    return (:head_loss, 1, if_alpha(alpha_m1, convex=true),
        f_alpha(alpha_m1, convex=true), df_alpha(alpha_m1, convex=true))
end

function head_loss_args(wm::AbstractMICPModel)
    alpha_m1 = get_alpha_min_1(wm)
    return (:head_loss, 1, f_alpha(alpha_m1, convex=true),
        df_alpha(alpha_m1, convex=true), d2f_alpha(alpha_m1, convex=true))
end

function head_loss_args(wm::AbstractNCNLPModel)
    alpha_m1 = get_alpha_min_1(wm)
    return (:head_loss, 1, f_alpha(alpha_m1, convex=false),
        df_alpha(alpha_m1, convex=false), d2f_alpha(alpha_m1, convex=false))
end

function primal_energy_args(wm::AbstractMICPModel)
    alpha = get_alpha_min_1(wm) + 1.0
    return (:primal_energy, 1, f_primal(alpha), df_primal(alpha), d2f_primal(alpha))
end

function dual_energy_args(wm::AbstractMICPModel)
    alpha = get_alpha_min_1(wm) + 1.0
    return (:dual_energy, 1, f_dual(alpha), df_dual(alpha), d2f_dual(alpha))
end

# By default, head loss is not defined by nonlinear registered functions.
function function_head_loss(wm::AbstractWaterModel) end

function function_head_loss(wm::AbstractNonlinearForms)
    JuMP.register(wm.model, head_loss_args(wm)...)
end

function function_head_loss(wm::MICPEWaterModel)
    JuMP.register(wm.model, head_loss_args(wm)...)
    JuMP.register(wm.model, primal_energy_args(wm)...)
    JuMP.register(wm.model, dual_energy_args(wm)...)
end
