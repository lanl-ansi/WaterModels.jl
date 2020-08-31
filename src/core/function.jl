#################################################################################
# This file defines the nonlinear head loss functions for water systems models. #
#################################################################################

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
            return x != 0.0 ? sign(x) * alpha * (1.0 + alpha) * (x*x)^(0.5*alpha - 0.5) : prevfloat(Inf)
        end
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

function head_loss_args(wm::CRDWaterModel)
    alpha_m1 = get_alpha_min_1(wm)
    return (:head_loss, 1, f_alpha(alpha_m1, convex=true),
        df_alpha(alpha_m1, convex=true), d2f_alpha(alpha_m1, convex=true))
end

function head_loss_args(wm::NCWaterModel)
    alpha_m1 = get_alpha_min_1(wm)
    return (:head_loss, 1, f_alpha(alpha_m1, convex=false),
        df_alpha(alpha_m1, convex=false), d2f_alpha(alpha_m1, convex=false))
end

# By default, head loss is not defined by nonlinear registered functions.
function function_head_loss(wm::AbstractWaterModel)
end

function function_head_loss(wm::AbstractNonlinearModel)
    JuMP.register(wm.model, head_loss_args(wm)...)
end
