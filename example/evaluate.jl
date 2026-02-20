using QEPOptimize
using BPGates

##
# Case A
##

# A case of no purification, verifying that the `PauliNoiseOp(i, f_in_to_pauli(network_fidelity)...)` works correctly

no_op_circuit = Individual([])

for network_fidelity in [0.7, 0.8, 0.9, 0.95]
    num_simulations = 1000000
    network_fidelity = 0.9
    noises = [NetworkFidelity(network_fidelity)]

    p1 = calculate_performance!(no_op_circuit; num_simulations, noises)

    @assert p1.success_probability == 1
    @assert abs(p1.error_probabilities[1] - network_fidelity) < 10/sqrt(num_simulations)
    @assert abs(p1.average_marginal_fidelity - network_fidelity) < 10/sqrt(num_simulations)

    p2 = calculate_performance!(no_op_circuit; num_simulations, noises, number_registers=3, purified_pairs=3)
    @assert p2.success_probability == 1
    @assert abs(p2.error_probabilities[1] - network_fidelity^3) < 10/sqrt(num_simulations)
    @assert abs(p2.average_marginal_fidelity - network_fidelity) < 10/sqrt(num_simulations)
end

##
# Case B
##

# The example circuit from figure 1 of https://quantum-journal.org/papers/q-2019-02-18-123/

network_fidelity = 0.9
noises = [NetworkFidelity(network_fidelity)]
number_registers = 2
purified_pairs = 1


# the CNOT gate is controlled on 1, targeting 2
# then we measure qubit 2
# many permutations of this circuit should give similar results
simple_purification = Individual([CNOTPerm(1,1,2,1), BellMeasure(1,2)])

num_simulations = 1000000

p = calculate_performance!(simple_purification; num_simulations, noises, number_registers, purified_pairs)

# check formula from appendix B
F = network_fidelity
q = (1-F)/3
prob = F^2+5q^2+2F*q
Fout = (F^2+q^2)/prob
@assert abs(p.purified_pairs_fidelity - Fout) < 10/sqrt(num_simulations)
@assert abs(p.success_probability - prob) < 10/sqrt(num_simulations)

##
# Case B'
##

# how we can do the same thing over a wide range of parameters

f_ins = [0.01; 0.05:0.05:0.95; 0.99; 0.999]
f_outs = Float64[]
probs = Float64[]
for f in f_ins
    noises = [NetworkFidelity(f)]
    simple_purification = Individual([CNOTPerm(1,1,2,1), BellMeasure(1,2)])
    p = calculate_performance!(simple_purification; num_simulations, noises, number_registers, purified_pairs)
    push!(f_outs, p.purified_pairs_fidelity)
    push!(probs, p.success_probability)
end

using CairoMakie
fig = Figure()
axF = Axis(fig[1,1])
axF.title = "F in vs F out"
lines!(axF, [0,1], [0,1], color=:gray50)
lines!(axF, f_ins, f_outs)
axP = Axis(fig[1,2])
axP.title = "F in vs P"
lines!(axP, f_ins, probs)
display(fig)
