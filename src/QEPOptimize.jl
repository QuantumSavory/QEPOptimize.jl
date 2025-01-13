module QEPOptimize

using BPGates
using BPGates: mctrajectory!, continue_stat, PauliNoise # TODO these should be exported by default

using Statistics: mean

using Random: randperm

export Individual, calculate_performance!, f_in_to_pauli, NetworkFidelity, NetworkPauliNoise, Population # TODO order these neatly

include("types.jl")
include("noises.jl") # TODO (low priority) this should be upstreamed to BPGates
include("performance_eval.jl")
include("mutate.jl")
include("evolve.jl")

end
