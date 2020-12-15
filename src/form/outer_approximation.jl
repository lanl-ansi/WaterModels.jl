function _get_owf_oa(q::JuMP.VariableRef, z::JuMP.VariableRef, q_hat::Float64, coeffs::Array{Float64})
    f = coeffs[1]*q_hat^3 + coeffs[2]*q_hat^2 + coeffs[3]*q_hat
    df = 3.0*coeffs[1]*q_hat^2 + 2.0*coeffs[2]*q_hat + coeffs[3]
    return f * z + df * (q - q_hat * z)
end


function _get_head_loss_oa(q::JuMP.VariableRef, q_hat::Float64, alpha::Float64)
    return q_hat^alpha + alpha * q_hat^(alpha - 1.0) * (q - q_hat)
end


function _get_head_loss_oa_binary(q::JuMP.VariableRef, z::Union{JuMP.VariableRef, JuMP.GenericAffExpr}, q_hat::Float64, exponent::Float64)
    return q_hat^exponent * z + exponent * q_hat^(exponent - 1.0) * (q - q_hat * z)
end


function _get_head_gain_oa(q::JuMP.VariableRef, z::JuMP.VariableRef, q_hat::Float64, curve_fun::Array{Float64})
    f = curve_fun[1]*q_hat^2 + curve_fun[2]*q_hat + curve_fun[3]
    df = 2.0 * curve_fun[1] * q_hat + curve_fun[2]
    return f * z + df * (q - q_hat * z)
end