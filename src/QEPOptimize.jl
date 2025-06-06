module QEPOptimize

using BPGates
using BPGates: mctrajectory!, continue_stat, PauliNoise # TODO these should be exported by default

using Makie

using Statistics: mean

using Random: randperm

using DataStructures: counter

using OhMyThreads: tmap,tmapreduce

export Individual, calculate_performance!, f_in_to_pauli, NetworkFidelity, NetworkPauliNoise, Population, # TODO order these neatly
    multiple_steps_with_history!,
    analyze_f_out_vs_f_in, plot_circuit_analysis, plot_fitness_history

include("types.jl")
include("noises.jl") # TODO (low priority) this should be upstreamed to BPGates
include("performance_eval.jl")
include("mutate.jl")
include("evolve.jl")
include("analysis.jl")

end
