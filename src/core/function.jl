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

function get_alpha_min_1(wm::AbstractWaterModel)
    alpha = [ref(wm, nw, :alpha) for nw in nw_ids(wm)]

    if !all(y -> y == alpha[1], alpha)
        Memento.error(_LOGGER, "Head loss exponents are different across the multinetwork.")
    else
        return alpha[1] - 1.0
    end
end

# By default, head loss is not defined by nonlinear registered functions.
function function_head_loss(wm::AbstractWaterModel) end

function function_head_loss(wm::AbstractCNLPModel)
    alpha = get_alpha_min_1(wm)
    f = JuMP.register(wm.model, :head_loss, 1, if_alpha(alpha, convex=true),
        f_alpha(alpha, convex=true), df_alpha(alpha, convex=true))
    #wm.fun[:head_loss] = (:head_loss, 1, if_alpha(alpha, convex=true),
    #    f_alpha(alpha, convex=true), df_alpha(alpha, convex=true))
end

function function_head_loss(wm::AbstractMICPModel)
    alpha = get_alpha_min_1(wm)
    f = JuMP.register(wm.model, :head_loss, 1, f_alpha(alpha, convex=true),
        df_alpha(alpha, convex=true), d2f_alpha(alpha, convex=true))
    #wm.fun[:head_loss] = (:head_loss, 1, f_alpha(alpha, convex=true),
    #    df_alpha(alpha, convex=true), d2f_alpha(alpha, convex=true))
end

function function_head_loss(wm::AbstractNCNLPModel)
    alpha = get_alpha_min_1(wm)
    f = JuMP.register(wm.model, :head_loss, 1, f_alpha(alpha, convex=false),
        df_alpha(alpha, convex=false), d2f_alpha(alpha, convex=false))
    #wm.fun[:head_loss] = (:head_loss, 1, f_alpha(alpha, convex=false),
    #    df_alpha(alpha, convex=false), d2f_alpha(alpha, convex=false))
end
