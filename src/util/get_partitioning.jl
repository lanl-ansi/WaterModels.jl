function get_partitioning(wm::GenericWaterModel,
                          resistance_indices::Dict{Int, Int},
                          solver::MathProgBase.AbstractMathProgSolver,
                          n::Int = wm.cnw)
    q, h = get_cvx_solution(wm, resistance_indices, solver, n)
    println(q)
    println(h)
end
