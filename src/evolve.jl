"""
Important: This calls sort and cull.

Execute one generation step of the genetic algorithm

Each generation, the population is optimized by applying mutation, selection, and crossover operations
while respecting constraints (such as gate connectivity and noise calibration)
"""
function step!(
    population::Population;
    max_ops::Int=5,
    number_registers::Int=2,
    purified_pairs::Int=1,
    num_simulations::Int=100,
    pop_size::Int=100,
    code_distance::Int=1,
    noises=[NetworkFidelity(0.9)],
    new_mutants::Int=10,
    p_drop=0.1,
    p_mutate=0.1,
    p_gain=0.1,
)
    # Mark existing individuals as survivors
    # Survivors ensure that some individuals are carried over unchanged, maintaining good solutions
    for indiv in population.individuals
        indiv.history = :survivor
    end

    # TODO parents and offspring circuits

    add_mutations!(population.individuals; max_ops, new_mutants, valid_pairs=1:number_registers) # TODO (low priority) decouple valid_pairs from number_registers

    # Sort the population by fitness and cull the excess individuals to maintain the population size
    simulate_and_sort!(
        population;
        num_simulations,
        purified_pairs,
        number_registers, # TODO (low priority) this should be by-default derived from `indiv`
        noises=[NetworkFidelity(0.9)] # TODO configurable noise
    )
    cull!(population, pop_size)
end

function multiple_steps_with_history!(
    population::Population, steps;
    step_callback=()->nothing,
    step_config... # same as step!
)
    fitness_history = Matrix{Float64}(undef, steps+1, step_config[:pop_size])
    fitness_history[1, :] = [i.fitness for i in population.individuals]
    transition_counts = []

    for i in 1:steps
        step!(population; step_config...)
        fitness_history[i+1,:] = [i.fitness for i in population.individuals]
        push!(transition_counts, counter([i.history for i in population.individuals]))
        step_callback()
    end

    transition_counts_keys = collect(union(keys.(transition_counts)...))
    transition_counts_matrix = Matrix{Int}(undef, length(transition_counts), length(transition_counts_keys))
    for (i, t) in enumerate(transition_counts)
        for (j, k) in enumerate(transition_counts_keys)
            transition_counts_matrix[i, j] = get(t, k, 0)
        end
    end
    return population, fitness_history, transition_counts_matrix, transition_counts_keys
end

"Adds mutated individuals to the given individual vector."
function add_mutations!(
    individuals::Vector{Individual};
    max_ops::Int=5,
    valid_pairs=1:1,
    new_mutants::Int=10,
    p_drop=0.1,
    p_mutate=0.1,
    p_gain=0.1,
    p_child=0.1,
    p_swap=0.1
)
    # a function that can be given to each thread, along with the required data, including an array to store the mutes
    function indiv_to_mutes(indiv)
        thread_mutants = Individual[]
        l = length(indiv.ops)

        for _ in 1:new_mutants
            # choose which mutations happen 
            _drop_op::Bool   = rand() < p_drop   && l > 0
            _gain_op::Bool   = rand() < p_gain   && l < max_ops
            _mutate::Bool = rand() < p_mutate && l > 0
            _child::Bool = length(individuals) > 1 && rand() < p_child && l > 0 
            _swap_op::Bool = rand() < p_swap && l > 1

            # do mutes 
            _drop_op && push!(thread_mutants, drop_op(indiv))
            _gain_op && push!(thread_mutants, gain_op(indiv; valid_pairs=valid_pairs))
            _mutate && push!(thread_mutants, mutate(indiv))
            _swap_op && push!(thread_mutants, swap_op(indiv.ops))

            if _child 
                # parents will be this indiv, and a random other one. Make sure it is not the same as this one 
                partner = rand(individuals)

                while partner === indiv
                    partner = rand(individuals)
                end
            
                partner_ops = partner.ops 

                # we have the partner now, make the child
                push!(thread_mutants, make_child(indiv.ops, partner_ops, max_ops))
            end
            
        end
        return thread_mutants
    end

    # Max threads will be set by the threads specified when running julia. ex) julia -t 16
    max_threads::Int = Threads.nthreads()

    # populate the thread_mutes by running the function on each thread
    mutants = tmapreduce(indiv_to_mutes,vcat,individuals; nchunks=max_threads) # TODO the reduce operation should be vcat

    ## add all mutes back to the individuals vector
    append!(individuals, mutants)
end

function sort_pop!(population::Population)
    # update the population with the sorted vector of individuals by fitness
    population.individuals = sort(population.individuals, by=indiv -> indiv.fitness, rev=true)
end

"Evaluate and Sort the individuals in descending order of fitness"
function simulate_and_sort!(
    population::Population;
    num_simulations::Int=100,
    purified_pairs::Int=1,
    number_registers::Int=2, # TODO (low priority) this should be by-default derived from `indiv`
    code_distance::Int=1,
    noises=[NetworkFidelity(0.9)]
)
    # calculate and update each individual's performance
    function update!(indiv)
        calculate_performance!(indiv;
            num_simulations,
            purified_pairs,
            number_registers,
            code_distance,
            noises)
        indiv.fitness = indiv.performance.purified_pairs_fidelity # TODO make it possible to select the type of fitness to evaluate
    end

    # Parallel processing for performance calculations. Max threads will be set by the threads specified when running julia. ex) julia -t 16
    max_threads::Int = Threads.nthreads()

    tmap(update!,population.individuals)
    
    sort_pop!(population)
end

"Reduce the population size to the target `population_size` (assumes pop is already sorted)"
function cull!(population::Population, population_size::Int)
    population.individuals = population.individuals[1:population_size]
end


"""
Initialize a population of quantum circuits, sort, and cull.
"""
function initialize_pop!(
    population::Population;
    start_ops::Int=10,
    start_pop_size::Int=1000,
    number_registers::Int=2,
    pop_size::Int=100,
    num_simulations::Int=100,
    purified_pairs::Int=1,
    code_distance::Int=1,
    noises=[NetworkFidelity(0.9)]
)
    valid_pairs=1:number_registers # TODO (low priority) decouple valid_pairs from number_registers

    number_registers >= 2 || throw(ArgumentError("number_registers must be >= 2"))

    for _ in 1:start_pop_size
        indiv = Individual(:random)
        for _ in 1:start_ops
            push!(indiv.ops, rand_op(valid_pairs))
        end
        push!(population.individuals, indiv)
    end

    simulate_and_sort!(
        population;
        num_simulations,
        purified_pairs,
        number_registers,
        code_distance,
        noises
    )
    cull!(population,pop_size)
end
