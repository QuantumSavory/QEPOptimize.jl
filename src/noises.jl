"""
    f_in_to_pauli(f_in)

Converts the `f_in` parameter to Pauli X, Y, and Z noise,
used in `calculate_performance!` to set up the initial noise.
"""
function f_in_to_pauli(f_in)
    px = py = pz = (1 - f_in) / 3
    return px, py, pz
end


"apply a given noise process to each operation in the circuit"
noisify_circuit(noise, circuit; number_registers) = noisify.((noise,), circuit)

"apply a given noise process to a given gate"
noisify(noise, operation) = operation # the default in absence of a defined method is that a given noise does not affect a given gate


"Pauli Noise experience during the establishment of raw Bell pairs over the network"
struct NetworkPauliNoise
    px::Float64
    py::Float64
    pz::Float64
end

"Measurement error: chance to report opposite result for coincidence measurement"
struct MeasurementError 
    p::Float64
end

"Thermal relaxation error to be applied per time step. Assuming same time step for all operations, same errors for all qubits"
struct T1T2Noise
    t1::Float64  # T1 time constant
    t2::Float64  # T2 time constant
    op_times::Dict{Type,Float64} # The time in seconds for each op. If not in the dict, time will be interpreted as zero.
    T1T2Noise(t1,t2,times) = new(t1,t2,times) 
end

const DEFAULT_OP_TIMES = Dict(
    BellMeasure => 1.0,
    CNOTPerm => 1.0,
    PauliNoiseBellGate => 1.0,
    NoisyBellMeasureNoisyReset => 1.0
)

T1T2Noise(t1::Float64,t2::Float64) = T1T2Noise(t1,t2,DEFAULT_OP_TIMES)
T1T2Noise(t1::Int,t2::Int) = T1T2Noise(t1,t2,DEFAULT_OP_TIMES)

# Retrive the time taken (seconds) during an operation, given the op_times dict, and op.
# This method is to get the top-level name of the type, (the wrapper) which makes it easy to setup a dict of op_times, as seen from the DEFAULT_OP_TIMES. Getting the wrapper effectively strips off the generic types like changing "PauliNoiseBellGate{CNOTPerm}" to "PauliNoiseBellGate" 
# In the future, if the generic types are specified in an op_times dict, this function will need to be changed to check "typeof(op)" instead of "typeof(op).name.wrapper"
function op_time(op_times::Dict{Type,Float64},op::Any) 
    try;    return op_times[typeof(op).name.wrapper] # Look for value in dict
    catch;  return 0.0 # (should be KeyError) no value found, so consider its time to be zero
    end
end
"A convenience constructor for unbiased Pauli network noise that results in a Bell pair of given fidelity"
NetworkFidelity(f) = NetworkPauliNoise(f_in_to_pauli(f)...)

function noisify_circuit(n::NetworkPauliNoise, circuit; number_registers)
    initial_noise = [PauliNoiseOp(i, n.px, n.py, n.pz) for i in 1:number_registers]
    return [initial_noise; noisify.((n,), circuit)]
end

noisify(n::NetworkPauliNoise, b::BellMeasure) = NoisyBellMeasureNoisyReset(b, 0, n.px, n.py, n.pz)
noisify(n::NetworkPauliNoise, b::NoisyBellMeasureNoisyReset) = NoisyBellMeasureNoisyReset(b.m, b.p, n.px, n.py, n.pz)

noisify(error::MeasurementError, b::NoisyBellMeasureNoisyReset) = NoisyBellMeasureNoisyReset(b.m, error.p, b.px, b.py, b.pz)
noisify(error::MeasurementError, b::BellMeasure) = NoisyBellMeasureNoisyReset(b, error.p, 0, 0, 0)

noisify(n::PauliNoise, c::CNOTPerm) = PauliNoiseBellGate(c, n.px, n.py, n.pz) # TODO reusing the QuantumClifford way of making noisy gates would be more consistent here

# "T1, T2 noise: find all 'idle' parts of the QC and fill them with T1 and T2 noise. No noise is applied during a op on that register, but only before and after until the end of the circuit
function noisify_circuit(error::T1T2Noise, circuit; number_registers)
    # tracks how much time has passed on each register since the last operation was applied, so if multiple registers on an incoming operation have the same time, they can receive no T1/T2 noise, but if they have different times, we can add T1/T2 noise to the one that has been idle, effectively simulating a parallel circuit operations with idle noise
    accounted_for_time = zeros(Float64, number_registers) 
    noisy_circuit = [] # to return
    @assert number_registers > 0 # this is assumed
    qubit_list = 1:number_registers # for use in indexing

    # Go through each operation
    for op in circuit 
        # get time that it would take 
        time_to_elapse = op_time(error.op_times, op)

        # make sure that the operations we care about are real actions on the qubit, not only noise
        # ie, do not add T1T2 noise between other T1T2 Noise! 
        # and ignore any other only-noise operations (ex. PauliNoiseOp)
        time_to_elapse == 0 && continue 

        @assert time_to_elapse > 0 # It should never be true that time is negative (Switch to unsigned?)

        ## --- Before applying operation ---

        # get the qubits it will act on
        acted_on_qubits = affectedqubits(op)

        # decide if we need to add padding before this operation 
        # if length(acted_on_qubits) == 1
        # one qubit is being acted on
        # we don't need to bring it up to speed with the other qubits yet, we only need to apply this operation on it. (Op applies in parallel)

        if length(acted_on_qubits) > 1
            # multiple qubits are being acted on
            # find the qubit with the maximum time accounted for, the other qubits must be 'brought up to speed' with it
            max_time = maximum([accounted_for_time[q] for q in 
            acted_on_qubits]) # get the maximum time accounted for

            # for every qubit that is not up to speed with the max time
            for q in acted_on_qubits
                if accounted_for_time[q] < max_time
                    # this qubit has not been accounted for up until now, so we must add T1 and T2 noise here now
                    
                    # get the error operations from time idling, t1, t2 values, and push them
                    idle_time = max_time - accounted_for_time[q]
                    λ₁, λ₂ = thermal_relaxation_error_rate(error.t1, error.t2, idle_time)
                    push!(noisy_circuit, T1NoiseOp(q,λ₁), T2NoiseOp(q,λ₂))

                    # now update its time
                    accounted_for_time[q] = max_time # set it to the maximum time, so that it is now in sync with the other qubits
                end
            end
        end

        # now all important qubits are up to speed, we can apply the operation 
        push!(noisy_circuit, op)

        # and update the time elapsed for this operation 
        for q in acted_on_qubits
            accounted_for_time[q] += time_to_elapse # update the time for all qubits that were acted on
        end

        ## --- After applying operation ---
    end

    # Now all operations have T1T2 noise padded before, and in-between them if needed. So now we need to bring all of the qubits up to speed with the end time.
    # Get the max time on a qubit (this is the time that it takes to execute the entire circuit)
    max_time = maximum(accounted_for_time)

    # Bring all qubits up to speed with the maximum time (pad with T1T2)
    for q in qubit_list
        if accounted_for_time[q] < max_time  
            # this qubit has not been accounted for up until now, so we must add T1 and T2 noise here now
            
            # get the error operations from time idling, t1, t2 values, and push them
            idle_time = max_time - accounted_for_time[q]
            λ₁, λ₂ = thermal_relaxation_error_rate(error.t1, error.t2, idle_time)
            push!(noisy_circuit, T1NoiseOp(q,λ₁), T2NoiseOp(q,λ₂))
        end
    end

    # TODO: make sure we can simplify/compactify_circuit with QuantumClifford

    return noisy_circuit
end
    
# thermal relaxation error from t1 and t2 values, and the time elapsed. From YipiaoWu/QuantumHardware
function thermal_relaxation_error_rate(t1, t2, gate_time) 
    λ₁ = 1 - exp(-gate_time/t1)
    t_ϕ = t1*t2 / (2*t1 - t2)
    λ₂ = 1 - exp(-gate_time/t_ϕ)
    return λ₁, λ₂
end


# These have been moved to BPGates.jl, and will be deleted after finishing the pull request open at https://github.com/QuantumSavory/BPGates.jl/pull/31
# if not using that version of BPGates, these will need to be uncommented.
# retrieve qubits involved in a gate
affectedqubits(gate::PauliNoiseBellGate) = (gate.g.idx1, gate.g.idx2)
affectedqubits(gate::NoisyBellMeasureNoisyReset) = (gate.m.sidx,)
affectedqubits(gate::BellMeasure) = (gate.sidx,)
affectedqubits(gate::CNOTPerm) = (gate.idx1, gate.idx2)
affectedqubits(gate)= ()  # Default case for gates that do not involve qubits

# TODO more types of noise should be implemented
    # gate_fidelity would turn CNOTPerm gates into gates wrapped into noise
    # T1/T2 noise will add noise that happens even during wait time
    # network_fidelity would turn BellMeasure into NoisyBellMeasureNoisyReset(...)
    # measurement_fidelity would do the same
    # Probably it would be best to define a `noisify(::AbstractNoise, ::AbstractOperation)::AbstractOperation`
    # Or even `noisify_circuit(::AbstractNoise, ::Vector{<:AbstractOperation})` so that we can also cover the T1, T2 noise
