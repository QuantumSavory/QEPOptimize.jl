using QEPOptimize
using QEPOptimize: initialize_pop!, step!, NetworkFidelity # TODO export these
using BPGates: PauliNoise # TODO re-export from QEPOptimize

using CairoMakie
using Quantikz: displaycircuit
using DataStructures: counter

##

config = (;
    num_simulations=1000, # needs to be large enough to resolve circuit noise but you can make smaller too # TODO anneal this and save old results
    number_registers=4, # do 3 for something faster
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

STEPS = 80
fitness_history = Matrix{Float64}(undef, STEPS+1, init_config.pop_size)
fitness_history[1, :] = [i.fitness for i in pop.individuals]
transition_counts = []

for i in 1:STEPS
    step!(pop; step_config...)
    fitness_history[i+1,:] = [i.fitness for i in pop.individuals]
    push!(transition_counts, counter([i.history for i in pop.individuals]))
end

transition_counts_keys = collect(union(keys.(transition_counts)...))
transition_counts_matrix = Matrix{Int}(undef, length(transition_counts), length(transition_counts_keys))
for (i, t) in enumerate(transition_counts)
    for (j, k) in enumerate(transition_counts_keys)
        transition_counts_matrix[i, j] = get(t, k, 0)
    end
end

##

fig = Figure(size=(600, 600))

# heatmap of fitness history
f_hist = fig[1,1]
a_hist = Axis(f_hist[1,1], xlabel="Generation", ylabel="Individual")
p_hist = plot!(a_hist, fitness_history[:,end:-1:begin], colorrange=(0.9, 1.0))
c_hist = Colorbar(f_hist[1,2], p_hist, label="Fitness")

# best and worst fitness over generations
f_bestworst = fig[2,1]
a_bestworst = Axis(f_bestworst[1,1], xlabel="Generation", ylabel="Fitness")
lines!(a_bestworst, maximum(fitness_history, dims=2)[:], label="Best")
lines!(a_bestworst, minimum(fitness_history, dims=2)[:], label="Worst")
axislegend(a_bestworst, position=:rb)

# type of circuit
f_type = fig[3,1]
a_type = Axis(f_type[1,1], xlabel="Generation", ylabel="Type")
for i in 1:size(transition_counts_matrix, 2)
    lines!(a_type, transition_counts_matrix[:,i], label=string(transition_counts_keys[i]))
end
axislegend(a_type, position=:rb)

fig

##

best_circuit = pop.individuals[1]
displaycircuit(best_circuit.ops)

##

f_ins = [0.01; 0.05:0.05:0.95; 0.99; 0.999]
f_outs = Float64[]
probs = Float64[]
num_simulations = 100000
for f in f_ins
    noises = [NetworkFidelity(f), PauliNoise(0.01/3, 0.01/3, 0.01/3)] # TODO plot multiple curves for different levels of Pauli Noise
    p = calculate_performance!(best_circuit; num_simulations, noises, config.number_registers, config.purified_pairs)
    push!(f_outs, p.purified_pairs_fidelity)
    push!(probs, p.success_probability)
end

fig = Figure()
axF = Axis(fig[1,1])
axF.title = "F in vs F out"
lines!(axF, [0,1], [0,1], color=:gray50)
lines!(axF, f_ins, f_outs)
axP = Axis(fig[1,2])
axP.title = "F in vs P"
lines!(axP, f_ins, probs)
display(fig)

##
