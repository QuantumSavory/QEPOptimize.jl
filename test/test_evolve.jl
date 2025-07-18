using TestItems

@testitem "simulation can run if pop_size changes"  begin
    using QEPOptimize:multiple_steps_with_history!,Population,initialize_pop!

    pop_size = 10
    pop = Population()
    initialize_pop!(pop; pop_size=pop_size)
    @test size(pop.individuals)[1] == pop_size

    multiple_steps_with_history!(pop,5;pop_size=pop_size)
    @test size(pop.individuals)[1] == pop_size

    pop_size = 5
    multiple_steps_with_history!(pop,5;pop_size=pop_size)
    @test size(pop.individuals)[1] == pop_size

    pop_size = 20
    multiple_steps_with_history!(pop,5;pop_size=pop_size)
    @test size(pop.individuals)[1] == pop_size
end

@testitem "Bad starting pop size still works" begin
    using QEPOptimize: initialize_pop!, multiple_steps_with_history!,Population
    # low new_mutants, low start_pop_size, high pop_size -> starting population is unideal (less than it should), but the simulation should still run without errors.
    config = (;
        num_simulations=1000, 
        number_registers=4,
        purified_pairs=1,
        pop_size = 100, # large
    )

    init_config = (;
        start_ops = 10,
        start_pop_size = 1, # small
        config...
    )

    step_config = (;
        max_ops = 15, 
        new_mutants = 5, # small
        p_drop = 0.1,
        p_mutate = 0.1,
        p_gain = 0.1,
        config...
    )

    pop = Population()
    initialize_pop!(pop; init_config...)
    _, fitness_history, transition_counts_matrix, transition_counts_keys = multiple_steps_with_history!(pop, 5; step_config...)

end