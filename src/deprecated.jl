# Function deprecation warnings.
# Can be removed in a breaking release after 09/01/2022.
function run_model(args...; kwargs...)
    @warn("`run_model` has been replaced with `solve_model`", maxlog = 1)
    solve_model(args...; kwargs...)
end


function run_des(args...; kwargs...)
    @warn("`run_des` has been replaced with `solve_des`", maxlog = 1)
    solve_des(args...; kwargs...)
end


function run_wf(args...; kwargs...)
    @warn("`run_wf` has been replaced with `solve_wf`", maxlog = 1)
    solve_wf(args...; kwargs...)
end


function run_mn_wf(args...; kwargs...)
    @warn("`run_mn_wf` has been replaced with `solve_mn_wf`", maxlog = 1)
    solve_mn_wf(args...; kwargs...)
end


function run_mn_wf_switching(args...; kwargs...)
    @warn("`run_mn_wf_switching` has been replaced with `solve_mn_wf_switching`", maxlog = 1)
    solve_mn_wf_switching(args...; kwargs...)
end


function run_owf(args...; kwargs...)
    @warn("`run_owf` has been replaced with `solve_owf`", maxlog = 1)
    solve_owf(args...; kwargs...)
end


function run_mn_owf(args...; kwargs...)
    @warn("`run_mn_owf` has been replaced with `solve_mn_owf`", maxlog = 1)
    solve_mn_owf(args...; kwargs...)
end