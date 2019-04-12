###############################################################################
# This file defines the nonlinear head loss functions for water systems models.
###############################################################################

function if_alpha(alpha::Float64)
    return function(x::Float64) return inv(2.0 + alpha) * (x*x)^(1.0 + 0.5*alpha) end
end

function f_alpha(alpha::Float64)
    # Stable version of x * abs(x)^(alpha)
    return function(x::Float64) return sign(x) * (x*x)^(0.5 + 0.5*alpha) end
end

function df_alpha(alpha::Float64)
    return function(x::Float64) return (1.0 + alpha) * (x*x)^(0.5*alpha) end
end

function d2f_alpha(alpha::Float64)
    return function(x::Float64) return sign(x) * alpha * (1.0 + alpha) * (x*x)^(0.5*alpha - 0.5) end
end

function function_f_alpha(wm::GenericWaterModel, alpha::Float64)
    JuMP.register(wm.model, :f_alpha, 1, f_alpha(alpha), df_alpha(alpha), d2f_alpha(alpha))
end

function function_if_alpha(wm::GenericWaterModel, alpha::Float64)
    JuMP.register(wm.model, :if_alpha, 1, if_alpha(alpha), f_alpha(alpha), df_alpha(alpha))
end
