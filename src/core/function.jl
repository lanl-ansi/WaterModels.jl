#################################################################################
# This file defines the nonlinear head loss functions for water systems models. #
#################################################################################


function _f_alpha(alpha::Float64; convex::Bool = false)
    if convex
        return function (x::Float64)
            return (x * x)^(0.5 + 0.5 * alpha)
        end
    else
        return function (x::Float64)
            return sign(x) * (x * x)^(0.5 + 0.5 * alpha)
        end
    end
end


function _df_alpha(alpha::Float64; convex::Bool = false)
    if convex
        return function (x::Float64)
            return (1.0 + alpha) * sign(x) * (x * x)^(0.5 * alpha)
        end
    else
        return function (x::Float64)
            return (1.0 + alpha) * (x * x)^(0.5 * alpha)
        end
    end
end


function _d2f_alpha(alpha::Float64; convex::Bool = false)
    if convex
        return function (x::Float64)
            return alpha * (1.0 + alpha) * (x * x)^(0.5 * alpha - 0.5)
        end
    else
        return function (x::Float64)
            return x != 0.0 ?
                   sign(x) * alpha * (1.0 + alpha) * (x * x)^(0.5 * alpha - 0.5) :
                   prevfloat(Inf)
        end
    end
end


function _get_alpha_min_1(wm::AbstractWaterModel)
    return _get_exponent_from_head_loss_form(wm.ref[:it][wm_it_sym][:head_loss]) - 1.0
end


function head_loss_args(wm::Union{NCDWaterModel,CRDWaterModel})
    alpha_m1 = _get_alpha_min_1(wm)
    return (
        :head_loss,
        1,
        _f_alpha(alpha_m1, convex = true),
        _df_alpha(alpha_m1, convex = true),
        _d2f_alpha(alpha_m1, convex = true),
    )
end


function head_loss_args(wm::NCWaterModel)
    alpha_m1 = _get_alpha_min_1(wm)
    
    return (
        :head_loss,
        1,
        _f_alpha(alpha_m1, convex = false),
        _df_alpha(alpha_m1, convex = false),
        _d2f_alpha(alpha_m1, convex = false),
    )
end


function _function_head_loss(wm::AbstractWaterModel)
    # By default, head loss is not defined by nonlinear registered functions.
end


function _function_head_loss(wm::AbstractNonlinearModel)
    JuMP.register(wm.model, head_loss_args(wm)...)
end
