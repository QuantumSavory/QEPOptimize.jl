
### Output to QASM notes
# Do operations on only one of the register in the pair 'alice only'
# Entanglement : custom op. 'generate ent.' Leave for the user to define
# Measurement : do measurement, and do nothing with result. After each 'non-final' measurement, bring entanglement back (do custom op again)
# Sparse gate, use custom op 'black box'
# Send msg for Sparsegates, stabilizer tab output


# TODO: testing with local QASM simulator

## Preparation of Bell pairs
# "create entanglement" default  
define_Φ⁺_qasm() = "// Create Psi plus Bell pair\ngate entangle a,b {\n\th a;\n\tcx a,b;\n}\n\n"
entangle_qasm(q1::Int,q2::Int) = "entangle q[$(q1)],q[$(q2)];\n" 
# Bell pair definitions left here for possible later use 
# ket:(|00⟩+|11⟩)/√2, Stabilizer:+XX +ZZ 
# Φ⁺_qasm(q1::Int,q2::Int) =  "h q[$(q1)];\ncx q[$(q1)],q[$(q2)];\n"

# # ket:(|00⟩-|11⟩)/√2, Stabilizer:-XX +ZZ 
# Φ⁻_qasm(q1::Int,q2::Int) = "x q[$(q1)];\n$(Φ⁺_qasm(q1,q2))"

# # ket:(|01⟩+|10⟩)/√2, Stabilizer:+XX -ZZ 
# Ψ⁺_qasm(q1::Int,q2::Int) = "x q[$(q2)];\n$(Φ⁺_qasm(q1,q2))"

# # ket:(|01⟩-|10⟩)/√2, Stabilizer:-XX -ZZ 
# Ψ⁻_qasm(q1::Int,q2::Int) = "$(Φ⁻_qasm(q1,q2))x q[$(q1)];\n"

# Default case
to_qasm(op) = "$(op) Not implemented\n"

# QuantumClifford Measurements 
# qasm: 'measure q -> b;' means measure qbit q in the z basis and record to bit b.
# TODO: Should each qbit return to the previous basis after being transformed? Ie, after performing hadamard and measure to measure in X, do inverse hadamard?
# Measurement in X basis. Hadamard, measure Z
to_qasm(op::sMX) = "h q[$(op.qubit)];\nmeasure q[$(op.qubit)] -> c[$(op.bit)];\n"
# Measurement in Y basis. S dagger, Hadamard, measure Z
to_qasm(op::sMY) = "sdg q[$(op.qubit)];\nh q[$(op.qubit)];\nmeasure q[$(op.qubit)] -> c[$(op.bit)];\n"
# Measurement in Z basis (normal measure)
to_qasm(op::sMZ) = "measure q[$(op.qubit)] -> c[$(op.bit)];\n"

""" 
Equivalent section from BPGates: 
const h = tHadamard
const p = tPhase
const hp = h*p
const ph = p*h
const good_perm_qc = ( # From the appendix of Optimized Entanglement Purification, but be careful with index notation being different
    (tId1,tId1), # TODO switch to symbolic gates
    (h*ph*ph,h*hp*hp*hp*hp),
    (h,h),
    (ph*ph,hp*hp*hp*hp),
    (ph,hp*hp),
    (h*ph,p*hp)
)"""
# translate this to qasm:
h(q::Int) = "h q[$q];\n"
p(q::Int) = "s q[$q];\n" 
hp(q::Int) = h(q) * p(q)
ph(q::Int) = p(q) * h(q)
id1(q::Int) = "id q[$q];\n"
"""Apply the two given functions in order, concat their outputs. Using the × unicode to avoid messing with Base module '*' functions.
So now we can say 
julia> (h×h)(1)
"h q[1];\nh q[1];\n"
"""
×(a,b) = (q::Int) -> a(q) * b(q)

const good_perm_qasm = (
    (id1,id1),
    (h×ph×ph,h×hp×hp×hp×hp),
    (h,h),
    (ph×ph,hp×hp×hp×hp),
    (ph,hp×hp),
    (h×ph,p×hp)
)

"""
    to_qasm(gate::CNOTPerm)
This is from BPGates:
    function toQCcircuit(gate::CNOTPerm)
    return [
        SparseGate(good_perm_qc[gate.single1][1], (gate.idx1*2-1,)),
        SparseGate(good_perm_qc[gate.single1][2], (gate.idx1*2,)),
        SparseGate(good_perm_qc[gate.single2][1], (gate.idx2*2-1,)),
        SparseGate(good_perm_qc[gate.single2][2], (gate.idx2*2,)),
        sCNOT(gate.idx1*2-1, gate.idx2*2-1),
        sCNOT(gate.idx1*2, gate.idx2*2)
    ]       
Instead of doing: BPGates.CNOTPerm -> SparseGate -> QASM,
and having to deal with sparsegate interpretation, we do this:
BPGates.CNOTPerm -> QASM, which means we will skip using toQCcircuit. 
We can hack this by making our own 'good_perm_qc' (good_perm_qasm) the same as how it is done in BPGates, except swapping out the gate definitions (h, p, hp, ph...) with the same QASM ones.
"""
function to_qasm(gate::CNOTPerm)
    return reduce(*,[
        good_perm_qasm[gate.single1][1](gate.idx1*2-1),
        good_perm_qasm[gate.single1][2](gate.idx1*2),
        good_perm_qasm[gate.single2][1](gate.idx2*2-1),
        good_perm_qasm[gate.single2][2](gate.idx2*2),
        to_qasm(sCNOT(gate.idx1*2-1, gate.idx2*2-1)),
        to_qasm(sCNOT(gate.idx1*2, gate.idx2*2))
    ])
end

"""
    to_qasm(gate::BellMeasure)
BellMeasure to qasm function, similar method to CNOTPerm
Equivalent function from BPGates to QuantumClifford:
function toQCcircuit(g::BellMeasure)
    meas = (sMX, sMY, sMZ)[g.midx]
    return [
        BellMeasurement([meas(g.sidx*2-1),meas(g.sidx*2)], g.midx==2)
        Reset(S"XX ZZ", [g.sidx*2-1,g.sidx*2])
    ]
end
Since we do not want to slow down quantum computation at all, (while waiting for classical computation) we will just perform a reset on/in the specified basis and qubit, and continue. We are not worrying about parity or comparing the results of the measurement.
"""
function to_qasm(gate::BellMeasure;re_entangle=nothing)
    meas = (sMX, sMY, sMZ)[gate.midx]

    qasm = reduce(*,[
        to_qasm(meas(gate.sidx*2-1)),
        to_qasm(meas(gate.sidx*2))
    ])
    if !isnothing(re_entangle) 
        qasm *= re_entangle(gate.sidx)
    end
    return qasm
end

# Normal CNOT
to_qasm(op::sCNOT) = "cx q[$(op.q1)], q[$(op.q2)]\n"

# To keep track of whether or not to re-entangle after some operations 
# should_re_entangle_after(op) = false
# should_re_entangle_after(op::BellMeasure) = true # only after measurement 
# others ?

"""
    to_qasm(circ, registers::Int,purified_pairs::Int;comments=true,entanglement=define_Φ⁺_qasm)
Convert Purification circuit to OPENQASM 3.0. Need to supply the circuit (vector/list of ops), total register amount, purified pair amount, and optionally ;comments=false to turn off comments.
"""
function to_qasm(circ, registers::Int,purified_pairs::Int;comments=true,entanglement=define_Φ⁺_qasm)
    # Should be some pure pairs, and need at least as many regs as pure pairs.
    @assert purified_pairs > 0 && purified_pairs <= registers

    # Get Header
    qasm_str = "OPENQASM 3.0;\ninclude \"stdgates.inc\";\n" # TODO: move to constant?

    # helper method to simply syntax: 'append to string'
    add(s::String) = qasm_str *= s

    # registers, and reset qubits
    add("qubit[$(registers*2)] q;\nreset q;\nbit[$(registers*2)] c;\n\n")


    # Define entanglement (defines the 'entangle' operation to be used to make pairs)
    add(entanglement())

    # Create bell pairs
    comments && add("// Creating Bell pairs\n")
    for i in 1:registers
        if comments
            # Mark purified pairs
            i <= purified_pairs ? add("// Purified pair $(i)\n") : add("// Non-purified pair $(i)\n")
        end
        add(entangle_qasm(2i-2, 2i-1)) # Create Bell pairs
    end

    # Add operations. 
    # TODO: Keep track of the operations, since we need to 're-entangle' physical qubits sometimes after measurement, but not at the end of the circuit.
    for (i,op) in enumerate(circ) # Now add the ops to the qasm output
        add("\n") # Add some space inbetween for readability
        comments && add("// Operation $(i): $(op)\n")
        # this is a hack, i need to fix this (issue around affectedqubits() )
        if isa(op,BellMeasure)
            # If it is a measurement, re-entangle (tell the function to manage )
            add(to_qasm(op;re_entangle=(pair)->entangle_qasm(2*pair-2, 2*pair-1)))
        else
            add(to_qasm(op))
        end
    end

    # pairs are now 'purified' so remind which ones. The logic results in one pair showing (0,1), and multiple pairs showing [(0,1),(2,3)...]
    add("\n// Entangled, purified pairs: $(purified_pairs == 1 ? (0, 1) : [(2i-2, 2i-1) for i in 1:purified_pairs])\n")

    comments && add("// Generated by QuantumSavory/QEPOptimize.jl\n")

    return qasm_str
end

# ignore any noise operations (which should not be in the circuit anyway, but might be) TODO: convert noisy ops to non-noisy ops
to_qasm(op::PauliNoise) = ""
