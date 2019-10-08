###############################################################################
# This file defines the nonlinear head loss functions for water systems models.
###############################################################################

function get_g_constants(alpha::Float64, delta::Float64)
    p::Float64 = alpha + 1.0

    a::Float64 = 15.0 * inv(8.0) * delta^(p - 1.0) +
        inv(8.0) * (p - 1.0) * p * delta^(p - 1.0) -
        7.0 * inv(8.0) * p * delta^(p - 1.0)

    b::Float64 = -5.0 * inv(4.0) * delta^(p - 3.0) -
        inv(4.0) * (p - 1.0) * p * delta^(p - 3.0) +
        5.0 * inv(4.0) * p * delta^(p - 3.0)

    c::Float64 = 3.0 * inv(8.0) * delta^(p - 5.0) +
        inv(8.0) * (p - 1.0) * p * delta^(p - 5.0) -
        3.0 * inv(8.0) * p * delta^(p - 5.0)

    return a, b, c
end

function g_alpha(alpha::Float64, delta::Float64)
    a, b, c = get_g_constants(alpha, delta)

    return function(x::Float64)
        return c*x*x*x*x*x + b*x*x*x + a*x
    end
end

function dg_alpha(alpha::Float64, delta::Float64)
    p::Float64 = alpha + 1.0
    a, b, c = get_g_constants(alpha, delta)

    return function(x::Float64)
        5.0*c * x*x*x*x + 3.0*b * x*x + a
    end
end

function d2g_alpha(alpha::Float64, delta::Float64)
    p::Float64 = alpha + 1.0
    a, b, c = get_g_constants(alpha, delta)

    return function(x::Float64)
        20.0*c * x*x*x + 6.0*b * x
    end
end

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

function function_head_loss(wm::AbstractWaterModel)
    # By default, head loss is not defined using nonlinear registered functions.
end

function function_head_loss(wm::AbstractCNLPModel)
    alpha = get_alpha_min_1(wm)
    f = JuMP.register(wm.model, :head_loss, 1, if_alpha(alpha, convex=true),
        f_alpha(alpha, convex=true), df_alpha(alpha, convex=true))
    wm.fun[:head_loss] = (:head_loss, 1, if_alpha(alpha, convex=true),
        f_alpha(alpha, convex=true), df_alpha(alpha, convex=true))
end

function function_head_loss(wm::AbstractMICPModel)
    alpha = get_alpha_min_1(wm)
    f = JuMP.register(wm.model, :head_loss, 1, f_alpha(alpha, convex=true),
        df_alpha(alpha, convex=true), d2f_alpha(alpha, convex=true))
    wm.fun[:head_loss] = (:head_loss, 1, f_alpha(alpha, convex=true),
        df_alpha(alpha, convex=true), d2f_alpha(alpha, convex=true))
end

function function_head_loss(wm::AbstractNCNLPModel)
    alpha = get_alpha_min_1(wm)
    f = JuMP.register(wm.model, :head_loss, 1, f_alpha(alpha, convex=false),
        df_alpha(alpha, convex=false), d2f_alpha(alpha, convex=false))
    wm.fun[:head_loss] = (:head_loss, 1, f_alpha(alpha, convex=false),
        df_alpha(alpha, convex=false), d2f_alpha(alpha, convex=false))
end
