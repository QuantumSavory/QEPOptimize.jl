# sweep_single.jl
# Run with:
# julia --project=<path to QEPOptimize env> sweep_single.jl
#
# Expects SLURM_ARRAY_TASK_ID in the environment.

using QEPOptimize
using QEPOptimize: initialize_pop!, multiple_steps_with_history!, NetworkFidelity
using BPGates: PauliNoise, T1Noise, T2Noise, MeasurementFlipNoise
using QuantumClifford: CircuitNoise
using Random
using Serialization

number_registers_range = [3, 4, 5]
network_fidelity_range = [0.85, 0.9, 0.95]
gate_error_range = [0.005, 0.01, 0.02]
measurement_flip_range = [0.0, 0.01, 0.02]

idle_noise_configs = [
    (; idle_type = :none, idle_mag = 0.0),
    (; idle_type = :t1, idle_mag = 0.01),
    (; idle_type = :t1, idle_mag = 0.03),
    (; idle_type = :t2, idle_mag = 0.01),
    (; idle_type = :t2, idle_mag = 0.03),
]

# Number of independent GA runs (different seeds) per parameter combination.
# Best-of-N gets picked later in an aggregation pass over results/.
const N_SEEDS = 3

param_grid = [
    (;
        number_registers = nr,
        network_fidelity = nf,
        gate_error = ge,
        measurement_flip = mf,
        idle_type = idle.idle_type,
        idle_mag = idle.idle_mag,
        seed_replicate = s,
    )
    for nr in number_registers_range
    for nf in network_fidelity_range
    for ge in gate_error_range
    for mf in measurement_flip_range
    for idle in idle_noise_configs
    for s in 1:N_SEEDS
]

const SIM = (;
    num_simulations = 1000,
    pop_size = 20,
    start_ops = 10,
    start_pop_size = 1000,
    max_ops = 15,
    new_mutants = 10,
    p_drop = 0.1,
    p_mutate = 0.1,
    p_gain = 0.1,
)

const N_STEPS = 80

task_id = parse(Int, ENV["SLURM_ARRAY_TASK_ID"])

if !(1 <= task_id <= length(param_grid))
    error(
        "SLURM_ARRAY_TASK_ID=$task_id is outside " *
        "the valid range 1:$(length(param_grid))"
    )
end

p = param_grid[task_id]

tag = "nr$(p.number_registers)" *
      "_nf$(p.network_fidelity)" *
      "_ge$(p.gate_error)" *
      "_mf$(p.measurement_flip)" *
      "_idle$(p.idle_type)$(p.idle_mag)" *
      "_seed$(p.seed_replicate)"

mkpath("results")
mkpath("results/failed")

println("[$task_id/$(length(param_grid))] $tag")
flush(stdout)

# Combine task_id with the seed replicate so different replicates of the
# same parameter combination actually get different random seeds, rather
# than task_id alone (which is already unique, so this is somewhat
# redundant, but makes the seed's dependence on the replicate explicit).
Random.seed!(hash((task_id, p.seed_replicate)))

try
    t0 = time()

    idle_noise_obj = if p.idle_type == :t1
        T1Noise(p.idle_mag)
    elseif p.idle_type == :t2
        T2Noise(p.idle_mag)
    else
        nothing
    end

    circuit_noise = CircuitNoise(
        single_qubit = PauliNoise(
            p.gate_error / 3,
            p.gate_error / 3,
            p.gate_error / 3,
        ),
        two_qubit = PauliNoise(
            p.gate_error / 3,
            p.gate_error / 3,
            p.gate_error / 3,
        ),
        idle_noise = idle_noise_obj,
        measurement = MeasurementFlipNoise(p.measurement_flip),
    )

    config = (;
        num_simulations = SIM.num_simulations,
        number_registers = p.number_registers,
        purified_pairs = 1,
        code_distance = 1,
        pop_size = SIM.pop_size,
        noises = [NetworkFidelity(p.network_fidelity)],
        circuit_noise,
    )

    init_config = (;
        start_ops = SIM.start_ops,
        start_pop_size = SIM.start_pop_size,
        config...,
    )

    step_config = (;
        max_ops = SIM.max_ops,
        new_mutants = SIM.new_mutants,
        p_drop = SIM.p_drop,
        p_mutate = SIM.p_mutate,
        p_gain = SIM.p_gain,
        config...,
    )

    pop = Population()

    initialize_pop!(
        pop;
        init_config...,
    )

    _, fitness_history, transition_counts_matrix, transition_counts_keys =
        multiple_steps_with_history!(
            pop,
            N_STEPS;
            step_config...,
        )

    elapsed = time() - t0
    best_circuit = pop.individuals[1]

    println(
        "[$task_id] finished in $(round(elapsed; digits=1)) s, " *
        "best fitness = $(round(best_circuit.fitness; digits=4))"
    )

    output_path = joinpath("results", "task$(task_id).jls")

    serialize(
        output_path,
        (;
            task_id,
            tag,
            params = p,
            fitness_history,
            transition_counts_matrix,
            transition_counts_keys,
            best_circuit,
            elapsed,
        ),
    )

    println("[$task_id] saved $output_path")

catch e
    bt = catch_backtrace()
    err_str = sprint(showerror, e, bt)

    println(stderr, "[$task_id] FAILED: $tag")
    println(stderr, err_str)
    flush(stderr)

    failure_path = joinpath("results", "failed", "task$(task_id).jls")
    serialize(
        failure_path,
        (;
            task_id,
            tag,
            params = p,
            error = err_str,
        ),
    )
    println("[$task_id] saved failure record to $failure_path")
    rethrow()
end