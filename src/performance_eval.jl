"""
Updates the individual's performance and returns it.

Simulates quantum operations using Monte Carlo trajectories to evaluate the performance of a given quantum purification circuit.
"""
function calculate_performance!(
    indiv::Individual;
    num_simulations::Int=100,
    purified_pairs::Int=1,
    number_registers::Int=1, # TODO (low priority) this should be by-default derived from `indiv`
    code_distance::Int=1,
    noises=[NetworkFidelity(0.9)]
)

    count_success = 0
    counts_marginals = zeros(Int,purified_pairs) # an array to find F₁, F₂, …, Fₖ (tracks how often each purified bell pair is in the desired state)
    counts_nb_errors = zeros(Int,purified_pairs+1) # an array to find P₀, P₁, …, Pₖ (tracks how often a given total number of errors happens) Careful with indexing it as it includes a P₀!

    noisy_purification_circuit = indiv.ops
    for n in noises
        noisy_purification_circuit = noisify_circuit(n, noisy_purification_circuit; number_registers)
    end

    for _ in 1:num_simulations
        purified_state, res = mctrajectory!(BellState(number_registers), noisy_purification_circuit)
        # If the circuit execution was reported as 'successful'
        if res == continue_stat
            count_success += 1
            err_count = 0
            for i in 1:purified_pairs # for each purified pair
                if purified_state.phases[2i-1] || purified_state.phases[2i] # checks whether an error has occurred based on binary representation in BPGates
                    err_count += 1
                else
                    counts_marginals[i] += 1 # tracks the i'th purified pair is in the desired state
                end
            end
            counts_nb_errors[err_count+1] += 1 # Careful with indexing it as it includes a P₀!
        end
    end

    if count_success == 0
        indiv.performance =  Performance(err_probs, 0,0, 0, 0)
        return indiv.performance
    end

    p_success = count_success    / num_simulations # proportion of successful simulations
    marginals = counts_marginals / count_success   # marginal fidelities of individual purified pairs
    err_probs = counts_nb_errors / count_success   # an array containing in each index i, how many errors occurred in i'th pair

    correctable_errors = div(code_distance  - 1, 2) # maximum number of correctable errors based on code distance after teleportation
    indiv_logical_qubit_fidelity = sum(err_probs[1:min(end, correctable_errors+1)]) # calculates the logical qubit fidelity by summing the probabilities of correctable errors

    indiv.performance =  Performance(err_probs, err_probs[1], indiv_logical_qubit_fidelity, mean(marginals), p_success)

    return indiv.performance
end
