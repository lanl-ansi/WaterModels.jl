function head_loss_integrated_hw_func(x)
    return 0.350631 * (x^2)^1.426
end

function head_loss_hw_func(x)
    return (x^2)^(0.926)
end

function head_loss_hw_prime(x)
    if x > 0.0
        return (1.852 * x) / (x^2)^(0.074)
    else
        return 0.0
    end
end

function head_loss_hw_prime_prime(x)
    if x > 0.0
        return -0.274096 / (x^2)^0.074 + 1.852 / (x^2)^0.074
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
