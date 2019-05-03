function find_initial_solution(wm::GenericWaterModel,
                               max_num_rounds::Int,
                               max_repair_iterations::Int,
                               fractional_decrease::Float64,
                               solver::MathProgBase.AbstractMathProgSolver,
                               n::Int = wm.cnw)
    # Initialize dictionaries to track resistance indices.
    connection_ids = collect(ids(wm, n, :connection))
    resistance_indices = Dict{Int, Int}(a => 1 for a in connection_ids) 
    best_resistance_indices = Dict{Int, Int}(a => 1 for a in connection_ids)

    # Set the initial solution such that arcs have minimum resistance.
    for (a, connection) in wm.ref[:nw][n][:connection]
        best_resistance_indices[a] = length(wm.ref[:nw][n][:resistance][a])
    end

    # Initialize the best objective value obtained thus far.
    best_objective_value = compute_objective(wm, best_resistance_indices, n)

    # Initialize the number of iterations used in the search procedure.
    num_rounds = 1

    while num_rounds <= max_num_rounds
        # Generate an initial solution by randomizing resistance selection.
        for (a, connection) in wm.ref[:nw][n][:connection]
            random_index = rand(1:length(wm.ref[:nw][n][:resistance][a]))
            resistance_indices[a] = random_index
        end

        # Compute the objective value corresponding to the solution above.
        objective_value = compute_objective(wm, resistance_indices, n)

        # While the objective value corresponding to the proposed solution has
        # not decreased by an amount that is "significant enough..."
        while objective_value > fractional_decrease * best_objective_value
            argmax_cost_reduction = nothing
            max_cost_reduction = 0.0

            # Find the local resistance change that will induce the largest
            # decrease in the current objective value.
            for (a, connection) in wm.ref[:nw][n][:connection]
                if resistance_indices[a] > 1
                    # Get the resistance index being used at this arc.
                    r = resistance_indices[a]

                    # Compute the corresponding reduction in cost.
                    cost_diff = wm.ref[:nw][n][:resistance_cost][a][r] -
                                wm.ref[:nw][n][:resistance_cost][a][r-1]
                    cost_reduction = connection["length"] * cost_diff

                    # If this is a better cost reduction, store it.
                    if cost_reduction > max_cost_reduction
                        max_cost_reduction = cost_reduction
                        argmax_cost_reduction = a
                    end
                end
            end

            if argmax_cost_reduction == nothing
                break
            end

            # Compute the resistance change and its corresponding objective.
            resistance_indices[argmax_cost_reduction] -= 1
            objective_value -= max_cost_reduction
        end

        # Repair the solution generated through the procedure above.
        repaired, resistance_indices = repair_solution(wm, resistance_indices,
                                                       max_repair_iterations,
                                                       best_objective_value,
                                                       solver, n)

        # Compute the objective value corresponding to the solution above.
        objective_value = compute_objective(wm, resistance_indices, n)

        # If the solution is feasible and the objective has decreased...
        if repaired && objective_value < best_objective_value
            best_objective_value = objective_value
            best_resistance_indices = deepcopy(resistance_indices)
        end

        # Increase the iteration counter.
        num_rounds += 1
    end

    # Return indices of the resistance choices.
    return best_resistance_indices
end
