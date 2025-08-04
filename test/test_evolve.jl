using TestItems

@testitem "Dynamic simulation count works in edge cases" begin
    using QEPOptimize
    using QEPOptimize: initialize_pop!, step!, NetworkFidelity # TODO export these
    using BPGates: PauliNoise # TODO re-export from QEPOptimize

    config = (;
        num_simulations=1000, 
        number_registers=4, 
        purified_pairs=1,
        code_distance=1,
        pop_size = 20,
        noises=[NetworkFidelity(0.9), PauliNoise(0.01/3, 0.01/3, 0.01/3)],
    )

    init_config = (;
        start_ops = 10,
        start_pop_size = 1000,
        config...
    )

    step_config = (;
        max_ops = 15,
        new_mutants = 10,
        p_drop = 0.1,
        p_mutate = 0.1,
        p_gain = 0.1,
        config...
    )

    # test that it still runs even though max simulation count might be less than sim count
    for max_s in [1,1500,3000]
        pop = Population()

        initialize_pop!(pop; init_config...)

        _, fitness_history, transition_counts_matrix, transition_counts_keys = multiple_steps_with_history!(pop, 20; max_simulations=max_s,step_config...);
    end

end
