# tune_hyperparams.jl
#
# Hyperparameter-tuning workflow for QEPOptimize.
#
# Run one phase at a time:
#
#   julia --project=. tune_hyperparams.jl
#
# Then, from the Julia REPL:
#
#   include("tune_hyperparams.jl")
#   tune_num_simulations()
#   tune_max_ops()
#   ...
#
# This file is not intended to be submitted as a full Slurm array blindly.
# Each candidate is repeated across several random seeds, and the saved
# results should be inspected before selecting the next set of values.

using QEPOptimize
using QEPOptimize:
    initialize_pop!,
    multiple_steps_with_history!,
    NetworkFidelity

using BPGates:
    PauliNoise,
    T1Noise,
    T2Noise,
    MeasurementFlipNoise

using QuantumClifford:
    CircuitNoise

using Random
using Serialization
using Statistics
using Sockets: gethostname


# --------------------------------------------------------------------------
# Fixed experiment decisions
# --------------------------------------------------------------------------

const BASELINE_METRIC = :logical_qubit_fidelity

const ALL_METRICS = [
    :logical_qubit_fidelity,
    :purified_pairs_fidelity,
    :average_marginal_fidelity,
]

# During optimizer tuning, keep structural parameters fixed so that
# GOOD/MODERATE/POOR differ only in physical noise conditions.
const GOOD = (
    number_registers = 4,
    network_fidelity = 0.95,
    gate_error = 0.005,
    measurement_flip = 0.0,
    idle_type = :none,
    idle_mag = 0.0,
    purified_pairs = 1,
)

const MODERATE = (
    number_registers = 4,
    network_fidelity = 0.90,
    gate_error = 0.01,
    measurement_flip = 0.01,
    idle_type = :t1,
    idle_mag = 0.01,
    purified_pairs = 1,
)

const POOR = (
    number_registers = 4,
    network_fidelity = 0.85,
    gate_error = 0.02,
    measurement_flip = 0.02,
    idle_type = :t2,
    idle_mag = 0.03,
    purified_pairs = 1,
)

const REPRESENTATIVE = (
    good = GOOD,
    moderate = MODERATE,
    poor = POOR,
)

# Held-out conditions. These may vary structural parameters because they are
# used only after hyperparameters have been frozen.
const UNSEEN_A = (
    number_registers = 3,
    network_fidelity = 0.92,
    gate_error = 0.015,
    measurement_flip = 0.005,
    idle_type = :t2,
    idle_mag = 0.02,
    purified_pairs = 1,
)

const UNSEEN_B = (
    number_registers = 5,
    network_fidelity = 0.87,
    gate_error = 0.008,
    measurement_flip = 0.02,
    idle_type = :t1,
    idle_mag = 0.02,
    purified_pairs = 3,
)

const UNSEEN = (
    a = UNSEEN_A,
    b = UNSEEN_B,
)

const RESULTS_DIR = joinpath(@__DIR__, "tuning_results")
mkpath(RESULTS_DIR)


# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------

function save_results(filename::AbstractString, value)
    path = joinpath(RESULTS_DIR, filename)
    serialize(path, value)
    println("Saved: $path")
    return path
end


"""Build a CircuitNoise instance from one condition."""
function make_circuit_noise(condition)
    idle_noise = if condition.idle_type == :t1
        T1Noise(condition.idle_mag)
    elseif condition.idle_type == :t2
        T2Noise(condition.idle_mag)
    elseif condition.idle_type == :none
        nothing
    else
        throw(ArgumentError(
            "Unsupported idle_type=$(condition.idle_type)"
        ))
    end

    return CircuitNoise(
        single_qubit = PauliNoise(
            condition.gate_error / 3,
            condition.gate_error / 3,
            condition.gate_error / 3,
        ),
        two_qubit = PauliNoise(
            condition.gate_error / 3,
            condition.gate_error / 3,
            condition.gate_error / 3,
        ),
        idle_noise = idle_noise,
        measurement = MeasurementFlipNoise(
            condition.measurement_flip
        ),
    )
end


"""Build the common QEPOptimize configuration."""
function make_config(
    condition;
    num_simulations::Int,
    pop_size::Int,
    code_distance::Int = 1,
)
    return (
        num_simulations = num_simulations,
        number_registers = condition.number_registers,
        purified_pairs = condition.purified_pairs,
        code_distance = code_distance,
        pop_size = pop_size,
        noises = [
            NetworkFidelity(condition.network_fidelity)
        ],
        circuit_noise = make_circuit_noise(condition),
    )
end


"""
Return the first generation after which all remaining fitness values stay
within `tol` of the final value.
"""
function plateau_index(
    history;
    tol::Float64 = 0.01,
)
    isempty(history) && return nothing

    target = history[end]
    threshold = tol * max(abs(target), 1e-9)

    for index in eachindex(history)
        stable = all(
            abs(value - target) <= threshold
            for value in history[index:end]
        )

        stable && return index
    end

    return nothing
end


"""
Extract the best individual without assuming population ordering.

Change `argmax` to `argmin` only if the selected metric is minimized.
"""
function best_individual(population)
    best_index = argmax(
        individual.fitness
        for individual in population.individuals
    )

    return population.individuals[best_index]
end


"""
Try to determine circuit length without hard-coding one representation.
Returns `missing` if no recognized field exists.
"""
function circuit_length(individual)
    if hasproperty(individual, :circuit)
        return length(getproperty(individual, :circuit))
    elseif hasproperty(individual, :operations)
        return length(getproperty(individual, :operations))
    elseif hasproperty(individual, :ops)
        return length(getproperty(individual, :ops))
    end

    return missing
end


"""Run one complete GA execution."""
function run_ga(
    condition;
    num_simulations::Int,
    max_performance_calcs::Int,
    max_ops::Int,
    pop_size::Int,
    new_mutants::Int,
    p_drop::Float64,
    p_mutate::Float64,
    p_gain::Float64,
    start_ops::Int,
    start_pop_size::Int,
    n_steps::Int,
    evolution_metric::Symbol,
    seed::Int,
    code_distance::Int = 1,
)
    start_ops <= max_ops || throw(ArgumentError(
        "start_ops=$start_ops must be <= max_ops=$max_ops"
    ))

    pop_size > 0 || throw(ArgumentError(
        "pop_size must be positive"
    ))

    new_mutants >= 0 || throw(ArgumentError(
        "new_mutants must be nonnegative"
    ))

    Random.seed!(seed)

    config = make_config(
        condition;
        num_simulations = num_simulations,
        pop_size = pop_size,
        code_distance = code_distance,
    )

    initialization_config = (
        start_ops = start_ops,
        start_pop_size = start_pop_size,
        config...,
    )

    step_config = (
        max_ops = max_ops,
        new_mutants = new_mutants,
        p_drop = p_drop,
        p_mutate = p_mutate,
        p_gain = p_gain,
        evolution_metric = evolution_metric,
        max_performance_calcs = max_performance_calcs,
        config...,
    )

    population = Population()

    initialize_pop!(
        population;
        initialization_config...,
    )

    start_time = time()

    _, fitness_history, transition_counts_matrix, transition_counts_keys =
        multiple_steps_with_history!(
            population,
            n_steps;
            step_config...,
        )

    elapsed = time() - start_time
    best = best_individual(population)

    return (
        seed = seed,
        condition = condition,
        best_fitness = best.fitness,
        best_circuit = best,
        best_circuit_length = circuit_length(best),
        fitness_history = fitness_history,
        plateau_generation = plateau_index(fitness_history),
        transition_counts_matrix = transition_counts_matrix,
        transition_counts_keys = transition_counts_keys,
        elapsed = elapsed,
        hostname = gethostname(),
        julia_threads = Threads.nthreads(),
        settings = (
            num_simulations = num_simulations,
            max_performance_calcs = max_performance_calcs,
            max_ops = max_ops,
            pop_size = pop_size,
            new_mutants = new_mutants,
            p_drop = p_drop,
            p_mutate = p_mutate,
            p_gain = p_gain,
            start_ops = start_ops,
            start_pop_size = start_pop_size,
            n_steps = n_steps,
            evolution_metric = evolution_metric,
            code_distance = code_distance,
        ),
    )
end


"""Run one candidate across several fixed seeds and summarize it."""
function run_repeated(
    condition;
    seeds = 1:3,
    kwargs...,
)
    runs = [
        run_ga(
            condition;
            kwargs...,
            seed = seed,
        )
        for seed in seeds
    ]

    fitnesses = getproperty.(runs, :best_fitness)
    elapsed_values = getproperty.(runs, :elapsed)

    plateau_values = [
        run.plateau_generation
        for run in runs
        if !isnothing(run.plateau_generation)
    ]

    circuit_lengths = [
        run.best_circuit_length
        for run in runs
        if !ismissing(run.best_circuit_length)
    ]

    return (
        runs = runs,
        seeds = collect(seeds),
        mean_fitness = mean(fitnesses),
        median_fitness = median(fitnesses),
        std_fitness = length(fitnesses) > 1 ?
            std(fitnesses) : 0.0,
        minimum_fitness = minimum(fitnesses),
        maximum_fitness = maximum(fitnesses),
        median_elapsed = median(elapsed_values),
        median_plateau_generation =
            isempty(plateau_values) ?
            nothing : median(plateau_values),
        median_circuit_length =
            isempty(circuit_lengths) ?
            missing : median(circuit_lengths),
    )
end


function print_summary(label, result)
    println(
        label,
        " median_fitness=",
        round(result.median_fitness; digits = 5),
        " mean=",
        round(result.mean_fitness; digits = 5),
        " std=",
        round(result.std_fitness; digits = 5),
        " median_elapsed=",
        round(result.median_elapsed; digits = 2),
        "s plateau=",
        result.median_plateau_generation,
        " circuit_length=",
        result.median_circuit_length,
    )
end


# --------------------------------------------------------------------------
# Step 1: full-run sensitivity to num_simulations
#
# Important: this does not isolate Monte Carlo estimator variance because
# every repetition performs a complete stochastic GA run. It measures how
# stable the final GA output is as num_simulations changes.
# --------------------------------------------------------------------------

function tune_num_simulations(
    condition = POOR;
    candidates = [100, 500, 1000, 2000, 5000],
    seeds = 1:5,
    max_performance_calcs = 10,
)
    println("=== STEP 1: num_simulations sensitivity ===")

    results = NamedTuple[]

    for num_simulations in candidates
        summary = run_repeated(
            condition;
            seeds = seeds,
            num_simulations = num_simulations,
            max_performance_calcs = max_performance_calcs,
            max_ops = 15,
            pop_size = 20,
            new_mutants = 10,
            p_drop = 0.1,
            p_mutate = 0.1,
            p_gain = 0.1,
            start_ops = 10,
            start_pop_size = 500,
            n_steps = 20,
            evolution_metric = BASELINE_METRIC,
        )

        print_summary(
            "num_simulations=$num_simulations",
            summary,
        )

        push!(
            results,
            (
                num_simulations = num_simulations,
                summary...,
            ),
        )
    end

    save_results(
        "step1_num_simulations.jls",
        results,
    )

    println(
        "Choose the smallest value where median fitness, variability, ",
        "and ranking behavior have stabilized."
    )

    return results
end


# --------------------------------------------------------------------------
# Step 1B: max_performance_calcs
# --------------------------------------------------------------------------

function tune_max_performance_calcs(
    condition = POOR;
    candidates = [5, 10, 25, 50],
    seeds = 1:3,
    num_simulations = 1000,
)
    println("=== STEP 1B: max_performance_calcs ===")

    results = NamedTuple[]

    for max_performance_calcs in candidates
        summary = run_repeated(
            condition;
            seeds = seeds,
            num_simulations = num_simulations,
            max_performance_calcs =
                max_performance_calcs,
            max_ops = 15,
            pop_size = 20,
            new_mutants = 10,
            p_drop = 0.1,
            p_mutate = 0.1,
            p_gain = 0.1,
            start_ops = 10,
            start_pop_size = 500,
            n_steps = 40,
            evolution_metric = BASELINE_METRIC,
        )

        print_summary(
            "max_performance_calcs=$max_performance_calcs",
            summary,
        )

        push!(
            results,
            (
                max_performance_calcs =
                    max_performance_calcs,
                summary...,
            ),
        )
    end

    save_results(
        "step1b_max_performance_calcs.jls",
        results,
    )

    return results
end


# --------------------------------------------------------------------------
# Step 2: max_ops
# --------------------------------------------------------------------------

function tune_max_ops(
    condition = POOR;
    candidates = [10, 15, 20, 25, 30],
    seeds = 1:3,
    num_simulations = 1000,
    pop_size = 20,
    n_steps = 80,
    max_performance_calcs = 10,
)
    println("=== STEP 2: max_ops ===")

    results = NamedTuple[]

    for max_ops in candidates
        start_ops = min(10, max_ops)

        summary = run_repeated(
            condition;
            seeds = seeds,
            num_simulations = num_simulations,
            max_performance_calcs =
                max_performance_calcs,
            max_ops = max_ops,
            pop_size = pop_size,
            new_mutants = 10,
            p_drop = 0.1,
            p_mutate = 0.1,
            p_gain = 0.1,
            start_ops = start_ops,
            start_pop_size = 1000,
            n_steps = n_steps,
            evolution_metric = BASELINE_METRIC,
        )

        print_summary(
            "max_ops=$max_ops",
            summary,
        )

        push!(
            results,
            (
                max_ops = max_ops,
                summary...,
            ),
        )
    end

    save_results(
        "step2_max_ops.jls",
        results,
    )

    println(
        "Choose the smallest max_ops where fitness plateaus and ",
        "best circuit lengths no longer repeatedly hit the ceiling."
    )

    return results
end


# --------------------------------------------------------------------------
# Step 3: population and generation budget
# --------------------------------------------------------------------------

function tune_population_budget(
    condition = POOR;
    candidates = [
        (pop_size = 20, new_mutants = 10, n_steps = 80),
        (pop_size = 20, new_mutants = 10, n_steps = 150),
        (pop_size = 40, new_mutants = 10, n_steps = 80),
        (pop_size = 40, new_mutants = 20, n_steps = 80),
        (pop_size = 80, new_mutants = 40, n_steps = 80),
    ],
    seeds = 1:3,
    num_simulations = 1000,
    max_ops = 15,
    max_performance_calcs = 10,
)
    println("=== STEP 3: population budget ===")

    results = NamedTuple[]

    for candidate in candidates
        summary = run_repeated(
            condition;
            seeds = seeds,
            num_simulations = num_simulations,
            max_performance_calcs =
                max_performance_calcs,
            max_ops = max_ops,
            pop_size = candidate.pop_size,
            new_mutants = candidate.new_mutants,
            p_drop = 0.1,
            p_mutate = 0.1,
            p_gain = 0.1,
            start_ops = min(10, max_ops),
            start_pop_size = 1000,
            n_steps = candidate.n_steps,
            evolution_metric = BASELINE_METRIC,
        )

        print_summary(
            "pop=$(candidate.pop_size) " *
            "mutants=$(candidate.new_mutants) " *
            "steps=$(candidate.n_steps)",
            summary,
        )

        push!(
            results,
            (
                candidate...,
                summary...,
            ),
        )
    end

    save_results(
        "step3_population_budget.jls",
        results,
    )

    return results
end


# --------------------------------------------------------------------------
# Step 4: mutation-probability combinations
# --------------------------------------------------------------------------

function tune_mutation_probs(
    condition = POOR;
    candidates = [
        (
            name = :balanced,
            p_drop = 0.10,
            p_mutate = 0.10,
            p_gain = 0.10,
        ),
        (
            name = :growth_biased,
            p_drop = 0.05,
            p_mutate = 0.10,
            p_gain = 0.20,
        ),
        (
            name = :shrink_biased,
            p_drop = 0.20,
            p_mutate = 0.10,
            p_gain = 0.05,
        ),
        (
            name = :mutation_biased,
            p_drop = 0.05,
            p_mutate = 0.20,
            p_gain = 0.05,
        ),
        (
            name = :aggressive,
            p_drop = 0.20,
            p_mutate = 0.20,
            p_gain = 0.20,
        ),
    ],
    seeds = 1:3,
    num_simulations = 1000,
    max_ops = 15,
    pop_size = 20,
    new_mutants = 10,
    n_steps = 80,
    max_performance_calcs = 10,
)
    println("=== STEP 4: mutation probabilities ===")

    results = NamedTuple[]

    for candidate in candidates
        summary = run_repeated(
            condition;
            seeds = seeds,
            num_simulations = num_simulations,
            max_performance_calcs =
                max_performance_calcs,
            max_ops = max_ops,
            pop_size = pop_size,
            new_mutants = new_mutants,
            p_drop = candidate.p_drop,
            p_mutate = candidate.p_mutate,
            p_gain = candidate.p_gain,
            start_ops = min(10, max_ops),
            start_pop_size = 1000,
            n_steps = n_steps,
            evolution_metric = BASELINE_METRIC,
        )

        print_summary(
            "mutation=$(candidate.name)",
            summary,
        )

        push!(
            results,
            (
                candidate...,
                summary...,
            ),
        )
    end

    save_results(
        "step4_mutation_probs.jls",
        results,
    )

    return results
end


# --------------------------------------------------------------------------
# Step 5: initial population
# --------------------------------------------------------------------------

function tune_initial_population(
    condition = POOR;
    start_ops_candidates = [5, 10, 15],
    start_pop_size_candidates = [250, 500, 1000, 2000],
    seeds = 1:3,
    num_simulations = 1000,
    max_ops = 15,
    pop_size = 20,
    new_mutants = 10,
    n_steps = 80,
    max_performance_calcs = 10,
)
    println("=== STEP 5: initial population ===")

    results = NamedTuple[]

    for start_ops in start_ops_candidates
        start_ops <= max_ops || continue

        for start_pop_size in start_pop_size_candidates
            summary = run_repeated(
                condition;
                seeds = seeds,
                num_simulations = num_simulations,
                max_performance_calcs =
                    max_performance_calcs,
                max_ops = max_ops,
                pop_size = pop_size,
                new_mutants = new_mutants,
                p_drop = 0.1,
                p_mutate = 0.1,
                p_gain = 0.1,
                start_ops = start_ops,
                start_pop_size = start_pop_size,
                n_steps = n_steps,
                evolution_metric = BASELINE_METRIC,
            )

            print_summary(
                "start_ops=$start_ops " *
                "start_pop=$start_pop_size",
                summary,
            )

            push!(
                results,
                (
                    start_ops = start_ops,
                    start_pop_size =
                        start_pop_size,
                    summary...,
                ),
            )
        end
    end

    save_results(
        "step5_initial_population.jls",
        results,
    )

    return results
end


# --------------------------------------------------------------------------
# Step 6: representative-condition cross-check
# --------------------------------------------------------------------------

function cross_check_conditions(
    frozen_settings;
    conditions = REPRESENTATIVE,
    metric = BASELINE_METRIC,
    seeds = 1:5,
)
    println("=== STEP 6: representative cross-check ===")

    results = NamedTuple[]

    for (condition_name, condition) in pairs(conditions)
        summary = run_repeated(
            condition;
            seeds = seeds,
            frozen_settings...,
            evolution_metric = metric,
        )

        print_summary(
            "condition=$condition_name",
            summary,
        )

        push!(
            results,
            (
                condition_name = condition_name,
                condition = condition,
                metric = metric,
                summary...,
            ),
        )
    end

    save_results(
        "step6_cross_check.jls",
        results,
    )

    return results
end


# --------------------------------------------------------------------------
# Step 7: unseen-condition validation
# --------------------------------------------------------------------------

function validate_unseen(
    frozen_settings;
    unseen = UNSEEN,
    metric = BASELINE_METRIC,
    seeds = 1:5,
)
    println("=== STEP 7: unseen validation ===")

    results = NamedTuple[]

    for (condition_name, condition) in pairs(unseen)
        summary = run_repeated(
            condition;
            seeds = seeds,
            frozen_settings...,
            evolution_metric = metric,
        )

        print_summary(
            "unseen=$condition_name",
            summary,
        )

        push!(
            results,
            (
                condition_name = condition_name,
                condition = condition,
                metric = metric,
                summary...,
            ),
        )
    end

    save_results(
        "step7_validate_unseen.jls",
        results,
    )

    return results
end


# --------------------------------------------------------------------------
# Step 8: metric cross-check
# --------------------------------------------------------------------------

function cross_check_metrics(
    frozen_settings;
    conditions = REPRESENTATIVE,
    metrics = ALL_METRICS,
    seeds = 1:3,
)
    println("=== STEP 8: metric cross-check ===")

    results = NamedTuple[]

    for metric in metrics
        for (condition_name, condition) in pairs(conditions)
            summary = run_repeated(
                condition;
                seeds = seeds,
                frozen_settings...,
                evolution_metric = metric,
            )

            print_summary(
                "metric=$metric " *
                "condition=$condition_name",
                summary,
            )

            push!(
                results,
                (
                    metric = metric,
                    condition_name =
                        condition_name,
                    condition = condition,
                    summary...,
                ),
            )
        end
    end

    save_results(
        "step8_cross_check_metrics.jls",
        results,
    )

    return results
end


# --------------------------------------------------------------------------
# Example usage
# --------------------------------------------------------------------------

# Run one phase at a time:
#
# tune_num_simulations()
# tune_max_performance_calcs()
# tune_max_ops()
# tune_population_budget()
# tune_mutation_probs()
# tune_initial_population()
#
# After selecting values:
#
# frozen = (
#     num_simulations = 1000,
#     max_performance_calcs = 10,
#     max_ops = 15,
#     pop_size = 20,
#     new_mutants = 10,
#     p_drop = 0.1,
#     p_mutate = 0.1,
#     p_gain = 0.1,
#     start_ops = 10,
#     start_pop_size = 1000,
#     n_steps = 80,
# )
#
# cross_check_conditions(frozen)
# validate_unseen(frozen)
# cross_check_metrics(frozen)
