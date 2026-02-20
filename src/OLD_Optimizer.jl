module Optimizer

# TODO:  add communication between qubits and hardware-specific errors from csv data to calculate_performance! Not the 'example error rates'

# TODO update exports, clean up what actually needs to be exported and what can remain internal to this file
# TODO: improve clarity of use and document how to use everything all together
# TODO: provide and document full usage of optimizer.jl within QEPO.jl
# TODO document Individual,population, and docs for all functions
# TODO: update abstractions and data types
# TODO: rewrite gain_op_with_constraints

"""
    Organize the noise in the ciruit. Some noise is in mutate, calculate performance, etc.

    Abstract away noise in initial state and in circuit

    Thermal relaxation notes for later:
    # add thermal relaxation noise T1 T2 into circuits
    # TODO: incorporate given hardware noise?
    # t1_avg, t2_avg, gate_times = 286e-6, 251e-6, 533e-9  # example for testing: average T1, T2 and t on 'ibmq_sherbrooke'
    # λ₁, λ₂ = thermal_relaxation_error_rate(t1_avg, t2_avg, gate_times)
    # noisy_ops = add_thermal_relaxation_noise(indiv.ops, λ₁, λ₂)

"""


"""
    Optimizer.jl - Genetic optimizer from YipiaoWu/QuantumHardware

    ***SETTING UP PARAMETERS
    In order to run the optimizer, these parameters are needed:

    1.  General optimizer parameters
        1a. Most important bits
        Stored in the 'Configuration' struct

        num_simulations::Int
        raw_bell_pairs::Int         # (n) incoming bell pairs
        purified_pairs::Int         # (k) outgoing bell pairs
        num_registers::Int          # amount of registers
        optimize_for::CostFunction  # optimization goal (see @enum CostFunction)
        max_gen::Int                # TODO more info
        max_ops::Int                # Limits the number of operations in each individual circuit


        1b.  Generic optimizer internals and quantum ciruit specifications
        Stored in the 'AdvancedConfiguration' struct

        code_distance::Int                              # for logical_qubit_fidelity
        communication_fidelity_in::Float64              # fidelity that being generated at communication qubits. 0.9 would be appropriate (f_in)
        population_size::Int                            #  target number of individuals in the population for each generation after initialization
        starting_pop_multiplier::Int                    #  A multiplier used to determine the size of the initial population
        starting_ops::Int                               # initial number of operations in each individual circuit of the population
        pairs::Int                                      # Number of parent pairs selected for breeding new individuals
        children_per_pair::Int                          # Number of offspring produced by each pair of parents
        mutants_per_individual_per_type::Int            # Number of mutations applied to each individual for each mutation type
        # Probabilities for various operations to be done during optimization:
        p_lose_operation::Float64
        p_add_operation::Float64
        p_swap_operation::Float64
        p_mutate_operation::Float64


        1c. Hardware specifications -> for IBM calibration data only, at the moment
        Stored in the 'HardwareConfiguration' struct

        calibration_data_path::String # Path to the csv file containing IBM calibration data
        calibration_data::Dict{Any,Any}
        valid_qubits::Array{Int}

    2. Formatting
        All of these parameters are stored like so, with the HardwareConfiguration and AdvancedConfiguration going inside of the Configuration struct, under the names
        hardware_config and advancd_config, respectively.

    ***RUNNING THE OPTIMIZER
    Once the config is set up and stored in a variable, create your population that will be used for the optimizer. For example this makes an empty population:
            population = Population()
    Now, you can call these functions:
            1. run_with_constraints!(population::Population, config::Configuration)
                - Runs the basic simulation and updates your population with the results
            2. run_with_constraints_history!(population::Population, config::Configuration)
                - Runs the simulation and updates your population with added data tracking, and plots the fidelities to see the optimizer's progress once finished running.

    The functions with ! are meant to change every part of the population, not the config.
"""

using Plots
using Base.Threads #  for multithreading
using Random
using Statistics
using LinearAlgebra
using BPGates # efficient representation of purification circuits
using QEPO.Configurable # import Configuration type from Configurable.jl
using QuantumClifford # general-purpose tools for Clifford circuits
using QuantumClifford.Experimental.NoisyCircuits
using BPGates: T1NoiseOp, T2NoiseOp
using OhMyThreads: tmap

# import Base: sort! # imports necessary for defining new methods of functions defined in Base

export generate_noisy_BellSwap_ops_for_individual, long_range_entanglement_generation!, Population, Performance, Individual, generate_valid_pairs,NoisyBellSwap, initialize_pop_with_constraints!,run_with_constraints_history!, ThreadData

### Genetic optimizer setup









# Create as struct for each thread's data used in calculate_performance
mutable struct ThreadData
    count_success::Int64
    counts_marginals::Vector{Int64}
    counts_nb_errors::Vector{Int64}
    err_count::Int64
end

### Quantum Circuit and qubit setup

"""convert two gate error rate to pauli noise representation"""
function p2_to_pauli(p2::Float64)
    px = py = pz = (1 - p2) / 4
    return px, py, pz
end


""" Convert p2 of two Bell pairs to Pauli representation """
function map_p2s_to_pauli(p2_A::Float64, p2_B::Float64)

    px_A = py_A = pz_A = p2_A/ 4        # pauli noise on bell pair A (two qubit)
    px_B = py_B = pz_B = p2_B/ 4        # pauli noise on bell pair B (two qubit)

    px = px_A + px_B - px_A * px_B
    py = py_A + py_B - py_A * py_B
    pz = pz_A + pz_B - pz_A * pz_B

    return px, py, pz
end


""" Convert pswap of two Bell pairs to Pauli representation """
function map_pswap_to_pauli(p2_A::Float64, p2_B::Float64)

    # calculate effective error probabilities of each swap gate
    # px_swap = 1 - (1-px)^3 -（1-px)*px^2 = 3*px - 4*px^2 + 2*px^3 ≈ 3*px
    # py_swap = 1 - (1-py)^3 -（1-py)*px^2 = 3*py - 4*py^2 + 2*py^3 ≈ 3*py
    # pz_swap = 1 - (1-pz)^3 -（1-pz)*px^2 = 3*pz - 4*pz^2 + 2*pz^3 ≈ 3*pz

    # infidities of swap gate is as 3 times as fidilities of cnot gate
    # For simplicity, higher-order terms can be ignored because p2 are small
    px_A = py_A = pz_A = p2_A/ 4 *3
    px_B = py_B = pz_B = p2_B/ 4 *3

    px = px_A + px_B - px_A * px_B
    py = py_A + py_B - py_A * py_B
    pz = pz_A + pz_B - pz_A * pz_B

    return px, py, pz
end

""" readout error probabilities of a Bell pair """
function pair_readout_error(eta_A::Float64, eta_B::Float64)

    eta = eta_A + eta_B - eta_A * eta_B    # eta = 1- (1-eta_A)(1-eta_B) = 1 - (1- eta_A - eta_B - eta_A * eta_B)

    return eta
end

""" Function to add thermal relaxation T1 and T2 noise to each gate in the circuit, including idle times """
# retrieve average t1 and t2 time
function avg_t1_t2(calibration_data, valid_qubits)
    t1_times = Float64[]
    t2_times = Float64[]
    for q in valid_qubits
        t1 = calibration_data[q][:T1]
        t2 = calibration_data[q][:T2]
        push!(t1_times, t1)
        push!(t2_times, t2)
    end
    t1_avg = mean(t1_times)
    t2_avg = mean(t2_times)
    return t1_avg, t2_avg
end

function thermal_relaxation_error_rate(t1, t2, gate_time) # experimental private function for internal use
    λ₁ = 1 - exp(-gate_time/t1)
    t_ϕ = t1*t2 / (2*t1 - t2)
    λ₂ = 1 - exp(-gate_time/t_ϕ)
    return λ₁, λ₂
end

# TODO this function does not add thermal noise to the "wait" times in between gates on the qubits that are not acted upon on the current timestep
function add_thermal_relaxation_noise(circuit, λ₁, λ₂) # experimental private function for intenral use
    max_steps = 100
    max_pairs = 4
    time_table = falses(max_steps, max_pairs)   # tracks the activity of each qubit over time
    thermal_noisy_circuit = []
    step = 1
    for gate in circuit
        # move to next time step if the current gate cannot be added at this step
        if !can_fill_table(time_table, gate, step)         # check if a gate can be added to the current time step
            add_idle_noise!(thermal_noisy_circuit, time_table, step, λ₁, λ₂)  # Add noise to idle qubits at each time step.
            step = step +1  # Move to the next step
        end
        push!(thermal_noisy_circuit, gate)                # Apply the original gate
        push_noise!(thermal_noisy_circuit, gate, λ₁, λ₂)  # Apply noise for active qubits in the gate
        fill_table!(time_table, gate, step)  # Mark qubits as active for this gate
    end

    # a workaround for slow performance -- this is a bad piece of code -- the proper way to fix this is to hook into the `compactify_circuit` capabilities of QuantumClifford or to just not use Julia ;)
    # UType = Union{PauliNoiseBellGate{CNOTPerm}, T1NoiseOp, PauliNoiseBellGate{BellSwap}, T2NoiseOp, NoisyBellMeasureNoisyReset}
    # UType = QuantumOperation
    # thermal_noisy_circuit = convert(Vector{UType}, thermal_noisy_circuit);
    return thermal_noisy_circuit
end

### helper functions
# Checks if a gate can be scheduled at the given step without conflicts
function can_fill_table(time_table, gate, step)
    qubits = get_qubits_involved(gate)
    for q in qubits
        if time_table[step, q]    # if is true, it means that qubit q is already engaged in another operation at this time step, so the gate cannot be scheduled
            return false  # Conflict: qubit is already active at this time step
        end
    end
    return true  # No conflicts; gate can be scheduled
end

# mark the qubits involved in the current gate as "active"
function fill_table!(time_table, gate, step)
    qubits = get_qubits_involved(gate)
    for q in qubits
        time_table[step, q] = true  # Mark qubit as active
    end
end

# retrieve qubits involved in a gate
get_qubits_involved(gate::PauliNoiseBellGate) = [gate.g.idx1, gate.g.idx2]
get_qubits_involved(gate::NoisyBellMeasureNoisyReset) = [gate.m.sidx]
get_qubits_involved(gate) = []

# Multimethod to apply noise for different gate types (in order to be extendable)
function push_noise!(circuit, gate::PauliNoiseBellGate{T}, λ₁, λ₂) where T
    push!(circuit, T1NoiseOp(gate.g.idx1, λ₁))
    push!(circuit, T2NoiseOp(gate.g.idx1, λ₂))
    push!(circuit, T1NoiseOp(gate.g.idx2, λ₁))
    push!(circuit, T2NoiseOp(gate.g.idx2, λ₂))
end

push_noise!(circuit, gate::NoisyBellMeasureNoisyReset, λ₁, λ₂) = []  # No thermal relaxation added to measurement
push_noise!(circuit, any::Any, λ₁, λ₂) = []
function add_idle_noise!(circuit, time_table, step, λ₁, λ₂)
    for q in 1:size(time_table, 2)
        if !time_table[step, q]  # If the qubit is idle at this time step
            push!(circuit, T1NoiseOp(q, λ₁))
            push!(circuit, T2NoiseOp(q, λ₂))
        end
    end
end


"""
    reset_population!(population,population_size::Int,starting_pop_multiplier::Int)

    Initialize individuals with empty operations and default performance
"""
function reset_population!(population,population_size::Int,starting_pop_multiplier::Int)
    reset_selection_history!(population)
    population.individuals=[
        Individual("random")
        for _ in 1:population_size * starting_pop_multiplier
    ]
end

### Setup algorithms


##################################################
###          initilize polulation              ###
##################################################

""" helper function to to map valid qubits to pairs """
function map_valid_qubits_to_pairs(valid_pair, valid_pairs)
    # create the dictionary where we are looking for the key
    valid_pair_to_num = Dict{Tuple{Int, Int}, Int}()
    for (i, pair) in enumerate(valid_pairs)
        valid_pair_to_num[pair] = i
        valid_pair_to_num[reverse(pair)] = i          # Ensures bidirectional access to pairs, treating (a, b) and (b, a) equivalently
    end
    return get(valid_pair_to_num, valid_pair, 0)       # retrieves the index of valid_pair from valid_pair_to_num, or 0 if valid_pair is not found
end

function map_num_to_valid_qubits(num::Int, valid_pairs)
    if num > 0 && num <= length(valid_pairs)
        return valid_pairs[num]
    end
    return (0, 0)
end


function generate_noisy_BellSwap_ops_for_individual(num_registers,valid_pairs,calibration_data)::Vector{Any}
    ##### Should num_gates:
    # num_gates = rand(1:get_starting_ops(advanced_config) - 1)  # Randomly determine the number of gates to include in each individual's circuit)
    ##### Be used here ??

    """ Create a sequence of BellSwap gate that will effectively move the qubits from the lowest index to the highest index in a structured manner"""
    swap_gates = []
    # for i in population.r:-1:2              # The outer loop ensures that for each qubit register from the second one to the topmost one

    for i in 2:1:num_registers
        for j in num_registers:-1:i          # The inner loop performs the swaps for the current qubit register i with all registers down to i-1
            push!(swap_gates, BellSwap(j, j-1))
        end
    end

    """ Wrap the noise of each BellSwap gate in a NoiseBellSwap based on the hardware's error rates """
    noisy_BellSwap = [
        begin
            pair_A, pair_B = map_num_to_valid_qubits(swap.idx1, valid_pairs), map_num_to_valid_qubits(swap.idx2, valid_pairs)
            # Are you indexing calibration_data with heterogenous keys? -- this is not a good idea in general
            # add two-qubit gate error rate from calibration data
            p2_A = calibration_data[pair_A].two_qubit_error
            p2_B = calibration_data[pair_B].two_qubit_error
            NoisyBellSwap(swap.idx1, swap.idx2, map_pswap_to_pauli(p2_A::Float64, p2_B::Float64)... )
        end
        for swap in swap_gates
    ]

    return noisy_BellSwap
end

""" Generate pairs in the order specified """
function generate_valid_pairs(valid_qubits)
    # Memoize this function so if called again for the same input, it will not need to recalculate
    memo = Dict{Array{Int, 1}, Array{Tuple{Int, Int}, 1}}()
    if haskey(memo, valid_qubits)
        return memo[valid_qubits]

    else
        valid_pairs = [(valid_qubits[i], valid_qubits[i + 1]) for i in 1:2:length(valid_qubits)-1]
        memo[valid_qubits] = valid_pairs
        return valid_pairs
    end
end

"""
    long_range_entanglement_generation!(population::Population,config::Configuration)

    For each individual in the population, set their operations to noisy BellSwaps.

    The sequence of BellSwap gates that will effectively move the qubits from the lowest index to the highest index in a structured manner, and wraps the gates in noisy bellswaps.
"""
function long_range_entanglement_generation!(population::Population,config::Configuration)
    reset_population!(population,config.advanced_config.population_size,config.advanced_config.starting_pop_multiplier,)

    """ Generate pairs in the order specified """
    valid_pairs = generate_valid_pairs(config.hardware_config.valid_qubits)

    """ Each individual in the population is processed in parallel to speed up the initialization process """
    Threads.@threads for indiv in population.individuals


        # num_gates = rand(1:get_starting_ops(advanced_config) - 1)  # Randomly determine the number of gates to include in each individual's circuit

        noisy_BellSwap = generate_noisy_BellSwap_ops_for_individual(config.num_registers, valid_pairs,config.hardware_config.calibration_data)

        indiv.ops = noisy_BellSwap

    end

end




"""
    reset_selection_history!(population::Population)

    Reset the selection history for the population
"""
function reset_selection_history!(population::Population)
    # TODO from Stefan: this list gets repeated frequently, probably it makes sense to put it in a "convenience" global variable (and as mentioned elsewhere, probably a list of symbols, not a list of strings)
    for hist in hist_list
        population.selection_history[hist] = Vector{Int64}()
    end

end


"""
    update_selection_history!(population::Population)

    Update the selection history for the individuals in the population
"""
function update_selection_history!(population::Population)
    for hist in hist_list
        push!(population.selection_history[hist], reduce(+, [1 for indiv in population.individuals if indiv.history  == hist], init=0)) # TODO from Stefan: `sum` is the same as `reduce(+)`
    end
end


##################################################
###          optimized process                 ###
##################################################

"""
    run_with_constraints!(population::Population, config::Configuration)

    Execution of a Genetic Algorithm Designed to Evolve a Population of Quantum Circuits
"""
function run_with_constraints!(population::Population, config::Configuration)
    # TODO from Stefan: the fact that there is `run` but also `run_with_constraints` kinda sounds like this can be written more neatly and simply if we use "multiple dispatch"

    initialize_pop_with_constraints!(population, config)

    # Sort and cull are now called automatically with the initialize_pop_with_constraints call.
    # simulate_and_sort!(population,config)  # Evaluate the performance and Sorts the individuals by fitness, placing the best-performing individuals at the top
    # cull!(population,config.advanced_config.population_size)  # Removes excess individuals to maintain the target population size

    for _ = 1:config.max_gen

        # Produce the next generation of individuals by performing selection, crossover, mutation, and other genetic operations
        step_with_constraints!(population,config.max_ops,config.hardware_config.valid_qubits,config.purified_pairs,config.num_registers,config.hardware_config.calibration_data,config.num_simulations,config.optimize_for,config.advanced_config)
        update_selection_history!(population)

        # Calculate performance for each individual in parallel using OhMyThreads.jl
        tmap(indiv -> calculate_performance!(indiv,
            config.num_simulations,
            config.purified_pairs,
            config.num_registers,
            config.optimize_for,
            config.advanced_config.code_distance,
            config.advanced_config.communication_fidelity_in,
            ), population.individuals)
    end

end

# alternative one when printing out the history
"""
    run_with_constraints_history!(population::Population, config::Configuration)

    Includes a plot for the history of the optimization for single F_in
TBW
"""
function run_with_constraints_history!(population::Population, config::Configuration)
    max_gen = config.max_gen

    max_purified_fidelities_gen_dic = zeros(max_gen,1)
    min_purified_fidelities_gen_dic = zeros(max_gen,1)
    ave_purified_fidelities_gen_dic = zeros(max_gen,1)

    initialize_pop_with_constraints!(population, config)
    # simulate_and_sort!(population,config)
    # cull!(population,population_size)

    for i = 1:config.max_gen
        # Mutate
        step_with_constraints!(population,config.max_ops,config.hardware_config.valid_qubits,config.purified_pairs,config.num_registers,config.hardware_config.calibration_data,config.num_simulations,config.optimize_for,config.advanced_config)
        update_selection_history!(population)

        # Calculate performance for each individual in parallel using OhMyThreads.jl
        performances = tmap(indiv -> calculate_performance!(indiv,
            config.num_simulations,
            config.purified_pairs,
            config.num_registers,
            config.optimize_for,
            config.advanced_config.code_distance,
            config.advanced_config.communication_fidelity_in
            ), population.individuals)

        purified_fidelities = [perf.purified_pairs_fidelity for perf in performances]
        # Filter out NaN values
        valid_fidelities = filter(!isnan, purified_fidelities)

        max_purified_fidelities_gen_dic[i]= maximum(valid_fidelities)
        min_purified_fidelities_gen_dic[i]= minimum(valid_fidelities)
        ave_purified_fidelities_gen_dic[i]= mean(valid_fidelities)

        # check running progress
        println("Running process: generation ",i)

    end

    # plot the history of the optimization for single F_in
    default(fontfamily="Times")
    generations = 1:max_gen
    plot(generations, max_purified_fidelities_gen_dic, label="Best purified pairs fidelity",lw=2,ylim=(0.7,1.01))
    plot!(generations, min_purified_fidelities_gen_dic, label="Worst purified pairs fidelity",lw=2,ylim=(0.7,1.01))
    plot!(generations, ave_purified_fidelities_gen_dic, label="Average purified pairs fidelity",lw=2,ylim=(0.7,1.01))
    xlabel!("Generation")
    ylabel!("Fitness")
    title!("Optimization History")

end


##########################################
##  (1) Parent Selection and Crossover   ##
##########################################

"""creates a new child individual from two parent individuals by combining their operations"""
function new_child(indiv::Individual, indiv2::Individual, max_ops::Int,num_registers::Int)::Individual
    if length(indiv2.ops) == 0
        return deepcopy(indiv) # No crossover if one of the parents has no operations
    elseif length(indiv.ops) == 0
        return deepcopy(indiv2)
    end
    new_indiv = deepcopy(indiv)

    ## filiter out the BellSwap gate
    # method 1
    ops1 = filter(op -> op isa PauliNoiseBellGate{CNOTPerm} || op isa NoisyBellMeasureNoisyReset, indiv.ops)
    ops2 = filter(op -> op isa PauliNoiseBellGate{CNOTPerm} || op isa NoisyBellMeasureNoisyReset, indiv2.ops)

    # method 2
    # i = indiv.r*(indiv.r-1)/2         # number of NoisyBellSwap gate is r*(r-1)/2
    # ops1 = indiv.ops[i+1:end]
    # ops1 = indiv.ops[i+1:end]

    # explore new circuit configurations
    # """With a 50% probability, the operations list of the first parent is reversed"""
    if rand() < 0.5
        ops1 = ops1[end:-1:1]
    end

    # """With a 50% probability, the operations list of the second parent is reversed"""
    if rand() < 0.5
        ops2 = ops2[end:-1:1]
    end

    # Potential problem: This selection may not be optimal because some valuable operations could be at the end of the parents' lists,
    # but the child might miss them due to the biased selection from the beginning only.

    # Randomly selects how many operations to take from each parent
    sample1 = 1:min(length(ops1), max_ops)
    if length(sample1) == 0
        return deepcopy(indiv2)
    end
    num_ops1 = rand(sample1)

    sample2 = 1:min(length(ops2), max_ops - num_ops1)
    if length(sample2) == 0
        return deepcopy(indiv)
    end
    num_ops2 = rand(sample2)

    # Combining the selected operations from both parents
    # TODO: fix types here?
    new_indiv.ops = vcat(
        # indiv.ops[1:num_registers*(num_registers-1)÷2],  # Include BellSwap gates
        ops1[1:num_ops1],
        ops2[1:num_ops2]
    )[1:min(end, max_ops)]  # Ensure the total number of operations does not exceed max_ops
    new_indiv.history =  "child"  # indicating that it was created through crossover

    return new_indiv

end


##########################################
##      (2) methods of mutations        ##
##########################################


"""
    swap_op_with_constraints(indiv::Individual)

    Swap two randomly selected operations, that can be cnot or measurement. Create new individual with these new ops and return it.
"""
function swap_op_with_constraints(indiv::Individual)::Individual
    if length(indiv.ops) < 2
        return indiv
    end
    new_indiv = deepcopy(indiv)
    ops = indiv.ops

    """ Randomly select two positions (cnot or measurement) """
    sample = [i for i in 1:length(ops) if isa(ops[i],PauliNoiseBellGate{CNOTPerm} ) || isa(ops[i], NoisyBellMeasureNoisyReset)]

    if length(sample) < 2
        return indiv
    end

    ind1 = rand(sample)

    # do not swap with the same operation
    sample2 = [i for i in 1:length(ops) if i != ind1 && (isa(ops[i],PauliNoiseBellGate{CNOTPerm} ) || isa(ops[i], NoisyBellMeasureNoisyReset))]

    if length(sample2) < 1
        return indiv
    end
    ind2 = rand(sample2)

    """ Get the operations at the selected positions """
    op1, op2 = ops[ind1], ops[ind2]

    """ Swap the operations """

    new_indiv.ops[ind1] = op2
    new_indiv.ops[ind2] = op1
    new_indiv.history = "swap_m"

    return new_indiv
end







#######################################################################################
##     renoise bell pair for different input fidilities when generating dataframe    ##
#######################################################################################

"""
    refresh_noise(n::PauliNoiseBellGate, f_in::Float64)

    reapply the noise parameters without considering p2
    No Change to Noise Parameters
TBW

"""
function refresh_noise(n::PauliNoiseBellGate, f_in::Float64)
    return PauliNoiseBellGate(n.g, n.px, n.py, n.pz)
end

"""
    refresh_noise(n::NoisyBellMeasureNoisyReset, f_in::Float64)

    The noise parameters are updated based on f_in
"""
function refresh_noise(n::NoisyBellMeasureNoisyReset, f_in::Float64)
    px, py, pz = f_in_to_pauli(f_in)
    return NoisyBellMeasureNoisyReset(n.m, n.p, px, py, pz)
end

"""
    refresh_noise(indiv::Individual, f_in::Float64)

    Reset and return an individual's performance and fitness, and refresh the noise of their operations.
"""
function refresh_noise(indiv::Individual, f_in::Float64)
    new_indiv = deepcopy(indiv)
    new_indiv.ops = [refresh_noise(op, f_in) for op in indiv.ops]
    new_indiv.performance = Performance([], 0, 0, 0, 0)
    new_indiv.fitness = 0.0
    return new_indiv
end

# TODO: get docs
"""
    total_raw_pairs(ops,num_registers,purified_pairs)

TBW
"""
function total_raw_pairs(ops,num_registers,purified_pairs)
    total = num_registers
    last_ops_reg = Set(1:num_registers)
    # TODO: fix the algorithm to not need to use deleition to count
    for op in reverse(ops)
        if isa(op, NoisyBellMeasureNoisyReset)
            t = op.m.sidx
            if t in last_ops_reg
                delete!(last_ops_reg, t)
                if t < purified_pairs
                    total += 1
                end
            else
                total += 1
            end
        elseif isa(op, PauliNoiseBellGate{CNOTPerm})
            for t in [op.g.idx1, op.g.idx2]
                delete!(last_ops_reg, t)
            end
        elseif isa(op, BellSwap)
            for t in [op.sidx1, op.sidx2]
                delete!(last_ops_reg, t)
            end
        end
    end

    return total
end

# TODO: check changed hash functionality
function Base.hash(indiv::Individual)
    return hash((Individual, [hash(op) for op in indiv.ops]))
end

function Base.hash(g::CNOTPerm)
    return hash((CNOTPerm, g.single1, g.single2, g.idx1, g.idx2))
end

function Base.hash(n::PauliNoiseBellGate{CNOTPerm})
    return hash((PauliNoiseBellGate{CNOTPerm}, n.px, n.py, n.pz))
end

function Base.hash(m::BellMeasure)
    return hash((BellMeasure, m.midx, m.sidx))
end

function Base.hash(n::NoisyBellMeasureNoisyReset)
    return hash((NoisyBellMeasureNoisyReset, hash(n.m), n.p, n.px, n.py, n.pz))
end


end
