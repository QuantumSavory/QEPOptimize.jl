using QEPOptimize
using QEPOptimize: initialize_pop!, step!, NetworkFidelity # TODO export these
using BPGates: PauliNoise # TODO re-export from QEPOptimize


##

config = (;
    num_simulations=1000, # needs to be large enough to resolve circuit noise but you can make smaller too # TODO anneal this and save old results
    number_registers=4, # do 3 for something faster
    purified_pairs=1,
    code_distance=1,
    pop_size = 20,
    noises=[NetworkFidelity(0.9), PauliNoise(0.01/3, 0.01/3, 0.01/3), 
            MeasurementError(0.1), T1T2Noise(1, 1)],
)

init_config = (;
    start_ops = 10,
    start_pop_size = 1000,
    config...
)

step_config = (;
    max_ops = 15, # do 5 for something faster
    new_mutants = 10,
    p_drop = 0.1,
    p_mutate = 0.1,
    p_gain = 0.1,
    config...
)

##

pop = Population()

initialize_pop!(pop; init_config...)

##

_, fitness_history, transition_counts_matrix, transition_counts_keys = multiple_steps_with_history!(pop, 80; step_config...);

##

# best_circuit = pop.individuals[1]

##
print("Finished")
