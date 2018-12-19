function head_loss_integrated_hw_func(x)
    return (250.0 * x * (x^2)^(463.0 / 500.0)) / 713.0
end

function head_loss_hw_func(x)
    return (x^2)^(0.926)
end

function head_loss_hw_prime(x)
    if x > 0.0
        return (463.0 * x) / (250.0 * (x^2)^(37.0 / 500.0))
    else
        return 0.0
    end
end

function head_loss_hw_prime_prime(x)
    if x > 0.0
        return 463.0 / (250.0 * (x^2)^(37.0 / 500.0)) -
               (17131.0*x^2) / (62500.0*(x^2)^(537.0 / 500.0))
    else
        return 0.0
    end
end

function function_head_loss_hw(wm::GenericWaterModel, n::Int = wm.cnw)
    JuMP.register(wm.model, :head_loss_hw, 1, head_loss_hw_func,
                  head_loss_hw_prime, head_loss_hw_prime_prime)
end

function function_head_loss_integrated_hw(wm::GenericWaterModel, n::Int = wm.cnw)
    JuMP.register(wm.model, :head_loss_integrated_hw, 1,
                  head_loss_integrated_hw_func,
                  head_loss_hw_func, head_loss_hw_prime)
end
