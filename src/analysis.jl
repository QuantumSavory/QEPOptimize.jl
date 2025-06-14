function analyze_f_out_vs_f_in(
    circuit::Individual;
    num_simulations::Int=100000,
    number_registers::Int=2,
    purified_pairs::Int=1,
    noises=[PauliNoise(0.01/3, 0.01/3, 0.01/3)],
    f_ins = [0.01; 0.05:0.05:0.95; 0.99; 0.999]
)
    f_outs = Float64[]
    probs = Float64[]
    for f in f_ins
        p = calculate_performance!(circuit; num_simulations, noises=[NetworkFidelity(f),noises...], number_registers, purified_pairs)
        push!(f_outs, p.purified_pairs_fidelity)
        push!(probs, p.success_probability)
    end
    return f_ins, f_outs, probs
end

function plot_circuit_analysis(
    circuit::Individual;
    num_simulations::Int=100000,
    number_registers::Int=2,
    purified_pairs::Int=1,
    noise_sets=[[PauliNoise(0.01/3, 0.01/3, 0.01/3)],[]],
    noise_set_labels=[join(string.(noises), " ") for noises in noise_sets],
    f_ins = [0.01; 0.05:0.05:0.95; 0.99; 0.999]
)
    fig = Figure()
    axF = Axis(fig[1,1])
    axF.title = "F in vs F out"
    lines!(axF, [0,1], [0,1], color=:gray50)
    axP = Axis(fig[1,2])
    axP.title = "F in vs P"
    for (noises, label) in zip(noise_sets, noise_set_labels)
        f_ins, f_outs, probs = analyze_f_out_vs_f_in(circuit; num_simulations, number_registers, purified_pairs, noises, f_ins)
        lines!(axF, f_ins, f_outs)
        lines!(axP, f_ins, probs; label)
    end
    axislegend(axP, "Noise Model", position=:lt)
    return fig
end

function plot_fitness_history(
    fitness_history, transition_counts_matrix, transition_counts_keys
)
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
    axislegend(a_type, position=:rb,labelsize=10,patchsize=(10,10))

    return fig
end

### Output to QASM 
# TODO: deal with SparseGates
# TODO: Implement BellMeasure, decide to "throw away" bad pairs or not. Not sure if this should be implemented, since QASM computation is split into sections of classical and quantum, so any classical checking will stop the quantum computation and could lead to wasting bell pairs.
    # // measure pair... (measure q[1] -> c[1])
    # // set up polarity of 'if' to be based on desired parity of BellMeasure 
    # // 'throw away operation, eg, reset and make a new pair
    # gate throw_away(index1,index2) {
    #   reset q[index1];
    #   reset q[index2];
    #   h q[index1];
    #   cx q[index1],q[index2];
    # }
    # if (c[1]==1) throw_away(0,1); // ?

# TODO: should Barriers be used? 
# TODO: testing with local QASM simulator

## Preperation of Bell pairs
# ket:(|00⟩+|11⟩)/√2, Stabilizer:+XX +ZZ 
Φ⁺_qasm(q1::Int,q2::Int) =  "h q[$(q1)];\ncx q[$(q1)],q[$(q2)];\n"

# ket:(|00⟩-|11⟩)/√2, Stabilizer:-XX +ZZ 
Φ⁻_qasm(q1::Int,q2::Int) = "x q[$(q1)];\n$(Φ⁺_qasm(q1,q2))"

# ket:(|01⟩+|10⟩)/√2, Stabilizer:+XX -ZZ 
Ψ⁺_qasm(q1::Int,q2::Int) = "x q[$(q2)];\n$(Φ⁺_qasm(q1,q2))"

# ket:(|01⟩-|10⟩)/√2, Stabilizer:-XX -ZZ 
Ψ⁻_qasm(q1::Int,q2::Int) = "$(Φ⁻_qasm(q1,q2))x q[$(q1)];\n"

# Default case
to_qasm(op) = "$(op) Not implemented\n"

# QuantumClifford Measurements 
# Measurement in X basis. Hadamard, measure Z
to_qasm(op::sMX) = "h q[$(op.qubit)];\nmeasure q[$(op.qubit)] -> c[$(op.bit)];\n"
# Measurement in Y basis. S dagger, Hadamard, measure Z
to_qasm(op::sMY) = "sdg q[$(op.qubit)];\nh q[$(op.qubit)];\nmeasure q[$(op.qubit)] -> c[$(op.bit)];\n"
# Measurement in Z basis (normal measure)
to_qasm(op::sMZ) = "measure q[$(op.qubit)] -> c[$(op.bit)];\n"

# CNOT permutation 
function to_qasm(op::CNOTPerm)
    # Convert BPGates.CNOTPerm into Vector{QuantumClifford.AbstractCliffordOperator}
    circ = toQCcircuit(op)
    # Then apply to_qasm() on each gate in vector, concat results
    return reduce(*,[to_qasm(QC_op) for QC_op in circ])
end

# SparseGates from CNOTPerm TODO
function to_qasm(op::SparseGate{Tableau{Vector{UInt8}, Matrix{UInt64}}})
    return "$(op) Not implemented\n"
end

# Normal CNOT
to_qasm(op::sCNOT) = "cx q[$(op.q1)], q[$(op.q2)]\n"

# Convert Purification circuit to OPENQASM 2.0. Need to supply the circuit (vector/list of ops), total register amount, purified pair amount, and optionally ;comments=false to turn off comments.
function to_qasm(circ, registers::Int,purified_pairs::Int;comments=true)
    # Get Header
    qasm_str = "OPENQASM 2.0;\ninclude \"qelib1.inc\";\n" # TODO: move to constant?

    # registers
    qasm_str *= "qreg q[$(registers*2)];creg c[$(registers*2)];\n"

    # Create bell pairs
    if comments
        qasm_str *= "// Creating Bell pairs\n"
    end
    for i in 1:registers
        if comments
            # Mark purified pairs
            if i <= purified_pairs
                qasm_str *= "// Purified pair $(i)\n"
            else
                qasm_str *= "// Non-purified pair $(i)\n"
            end
        end
        qasm_str *= Φ⁺_qasm(2i-2, 2i-1) # Create Bell pairs
    end

    # Add operations
    for (i,op) in enumerate(circ)
        if comments 
            qasm_str *= "// Operation $(i): $(op)\n"
        end
        qasm_str *= to_qasm(op)
    end

    if comments
        qasm_str *= "Genetated by QuantumSavory/QEPOptimize.jl\n"
    end
    return qasm_str
end

# ignore any noise operations (which should not be in the circuit anyway, but might be) TODO: convert noisy ops to non-noisy ops
to_qasm(op::PauliNoise) = ""

### Convert to Stabilizer view TODO 
# Very rough beginning draft of this function 
function to_stabilizer(circ,registers)
    state = BellState(registers)
    stab = MixedDestabilizer(state)
    output = string(stab)
    for op in circ
        # Get Clifford version of op
        for cliff_op in toQCcircuit(op)
                # Apply op to stabilizer
                if isa(cliff_op,BellMeasurement)
                    for m_op in cliff_op.measurements
                        apply!(stab,m_op)
                    end
                # TODO add more ops here/switch to multiple dispatch on to_stabilizer() calls
                # ie, define to_stabilizer function for any possible op, then combine results of calls
                else
                    apply!(stab,cliff_op)
                end
                # Append output
                output *= "\n$(cliff_op)\n$(string(stab))"
        end
    end
    return output
end
