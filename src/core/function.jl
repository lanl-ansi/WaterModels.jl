

# This file defines the nonlinear head loss functions for water systems models. #
#################################################################################

function _get_alpha(wm::AbstractWaterModel)
    head_loss_form = wm.ref[:it][wm_it_sym][:head_loss]
    alphas = uppercase(head_loss_form) == "H-W" ? 1.852 : 2.0
    # alphas = [ref(wm, nw, :alpha) for nw in nw_ids(wm)]
    if any(x -> x != alphas[1], alphas)
        Memento.error(_LOGGER, "Head loss exponents are different across the multinetwork.")
    end
    return alphas[1]
end

function head_loss(wm::Union{NCDWaterModel,CRDWaterModel}, x)
    p = _get_alpha(wm)
    # An expression equivalent to abspower(x, p) = abs(x)^p
    # return JuMP.@expression(wm.model, ifelse(x > 0, x^p, (-x)^p))
    return JuMP.@expression(wm.model, (x^2)^(p / 2))
end

function head_loss(wm::AbstractNonlinearModel, x)
    p = _get_alpha(wm)
    # An expression equivalent to signpower(x, p) = sign(x) * abs(x)^p
    return JuMP.@expression(wm.model, ifelse(x > 0, x^p, -(-x)^p))
end
