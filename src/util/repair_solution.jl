function repair_solution(wm::GenericWaterModel,
                         resistance_indices::Dict{Int, Int},
                         max_iterations::Int,
                         best_objective_value::Float64,
                         solver::MathProgBase.AbstractMathProgSolver,
                         n::Int = wm.cnw)
    # Initialize parameters used within the while loop.
    num_iterations = 1 # Number of iterations used throughout the repair loop.
    progress = true # Indicates whether changes were made to resistances.
    repaired = false # Indicates whether resistance choices have been repaired.

    while num_iterations <= max_iterations && !repaired && progress
        # Get the solution to CVXNLP for the current set of resistances.
        q, h = get_cvx_solution(wm, resistance_indices, solver)

        # Get bound test results and determine if the solution is feasible.
        qlb, qub, hlb, hub = check_solution_bounds(wm, q, h, resistance_indices, n)
        solution_is_feasible = all([all(values(qlb)), all(values(qub)),
                                    all(values(hlb)), all(values(hub))])

        if solution_is_feasible
            repaired = true # Solution does not need to be repaired.
        else
            progress = false # Repair progress has not yet been made.

            # Repair resistances where maximum flow violations have occurred.
            for (a, connection) in wm.ref[:nw][n][:connection]
                # Total number of possible resistances along arc a.
                n_r = length(wm.ref[:nw][n][:resistance][a])

                # If a maximum flow violation has occurred along this arc...
                if !qub[a] && resistance_indices[a] < n_r
                    resistance_indices[a] += 1
                    progress = true
                end
            end

            if !progress # If progress has not yet been made...
                # Repair resistances where minimum head has been violated.
                for (a, connection) in wm.ref[:nw][n][:connection]
                    # Get the total number of possible resistances.
                    n_r = length(wm.ref[:nw][n][:resistance][a])

                    if resistance_indices[a] < n_r
                        # Get indices for the two nodes connecting the arc.
                        i, j = get_node_ids(connection)

                        # Compute if a downstream node violation has occurred.
                        violation = (q[a] >= 0.0 && !hlb[j] && hlb[i]) ||
                                    (q[a] < 0.0 && !hlb[i] && hlb[j])

                        # If a downstream node violation has occurred...
                        if violation && resistance_indices[a] < n_r
                            resistance_indices[a] += 1
                            progress = true
                        end
                    end
                end
            end
        end

        # Compute the objective value for the selected resistance indices.
        objective_value = compute_objective(wm, resistance_indices, n)

        if objective_value >= best_objective_value
            # If the objective has not improved, progress has not been made.
            progress = repaired = false
        end

        # Increase the iteration counter.
        num_iterations += 1
    end

    # Return both the status of whether or not the solution has been repaired
    # (true or false) as well as the indices of the resistance choices.
    return repaired, resistance_indices
end
