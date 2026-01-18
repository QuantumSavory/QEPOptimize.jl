using TestItems

@testitem "performance averaging" begin
    using QEPOptimize
    using QEPOptimize:Performance

    low = Performance([0,0,1],0,0,0,0,1)
    high = Performance([1,0,0],1,1,1,1,1)
    mid = low + high
    @test mid == Performance([0.5,0,0.5],0.5,0.5,0.5,0.5,2)

    # works with +=
    low = Performance([0,0,1],0,0,0,0,1)
    high = Performance([1,0,0],1,1,1,1,1)

    low += high
    @test low == mid
end

@testitem "Performance averages deal with low sim count/circuit failures" begin 
    using QEPOptimize
    using QEPOptimize: initialize_pop!, step!, NetworkFidelity 
    using BPGates: PauliNoise 

    config = (;
        num_simulations=5, # needs to be large enough to resolve circuit noise but you can make smaller too # TODO anneal this and save old results
        number_registers=4, # do 3 for something faster
        purified_pairs=1,
        code_distance=1,
        pop_size = 20,
        noises=[NetworkFidelity(0.5), PauliNoise(0.01/3, 0.01/3, 0.01/3)],
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

    pop = Population()

    initialize_pop!(pop; init_config...)

    _, fitness_history, transition_counts_matrix, transition_counts_keys = multiple_steps_with_history!(pop, 80; step_config...)

end

@testitem "performance averages deal with bad error_probs" begin
    using QEPOptimize
    using QEPOptimize:Performance

    # defaults to the second performance if error probs are different
    normal = Performance([0,0,1],0,0,0,0,1)
    new_err_probs = Performance([1,0,0,0],1,1,1,1,1)
    result = normal + new_err_probs
    @test result == Performance([1,0,0,0],1,1,1,1,1)

    # defaults to the second err_probs if no error probs, still averages
    normal = Performance([],0,0,0,0,1)
    new_err_probs = Performance([1,0,0,0],1,1,1,1,1)
    result = normal + new_err_probs
    @test result == Performance([1,0,0,0],0.5,0.5,.5,.5,2)
    
end
