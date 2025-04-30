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
)
    mutants = Individual[]

    # For every mutation per individual, up to the limit
    for old_indiv in individuals
        for _ in 1:new_mutants
            indiv = copy(old_indiv)
            l = length(indiv.ops)
            _drop_op::Bool   = rand() < p_drop   && l > 0
            _gain_op::Bool   = rand() < p_gain   && l < max_ops
            # swap_op::Bool   = TODO
            _mutate::Bool = rand() < p_mutate && l > 0

            _drop_op && push!(mutants, drop_op(indiv))
            _gain_op && push!(mutants, gain_op(indiv; valid_pairs))
            _mutate && push!(mutants, mutate(indiv))
        end
    end

    ## add all children back to the individuals vector
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
    #=Threads.@threads=# for indiv in population.individuals
        calculate_performance!(indiv;
            num_simulations,
            purified_pairs,
            number_registers,
            code_distance,
            noises)
        indiv.fitness = indiv.performance.purified_pairs_fidelity # TODO make it possible to select the type of fitness to evaluate
    end
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
