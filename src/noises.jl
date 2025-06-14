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
end

"A convenience constructor for unbiased Pauli network noise that results in a Bell pair of given fidelity"
NetworkFidelity(f) = NetworkPauliNoise(f_in_to_pauli(f)...)

function noisify_circuit(n::NetworkPauliNoise, circuit; number_registers)
    initial_noise = [PauliNoiseOp(i, n.px, n.py, n.pz) for i in 1:number_registers]
    return [initial_noise; noisify.((n,), circuit)]
end

noisify(n::NetworkPauliNoise, b::BellMeasure) = NoisyBellMeasureNoisyReset(b, 0, n.px, n.py, n.pz)
noisify(n::NetworkPauliNoise, b::NoisyBellMeasureNoisyReset) = NoisyBellMeasureNoisyReset(b.m, b.p, n.px, n.py, n.pz)

function noisify_circuit(error::MeasurementError, circuit; number_registers)
    return map(op -> noisify(error, op), circuit)
end

noisify(error::MeasurementError, b::NoisyBellMeasureNoisyReset) = NoisyBellMeasureNoisyReset(b.m, error.p, b.px, b.py, b.pz)
noisify(error::MeasurementError, b::BellMeasure) = NoisyBellMeasureNoisyReset(b, error.p, 0, 0, 0)

noisify(n::PauliNoise, c::CNOTPerm) = PauliNoiseBellGate(c, n.px, n.py, n.pz) # TODO reusing the QuantumClifford way of making noisy gates would be more consistent here

# These will provide an estimated amount of time each operation will take, for use in `noisify_circuit(error::T1T2Noise, circuit::Vector{Any}; number_registers)` 
# TODO: units/actual data
op_time(op::CNOTPerm) = 1.0  
op_time(op::BellMeasure) = 1.0  
op_time(op::NoisyBellMeasureNoisyReset) = 1.0  
op_time(op::PauliNoiseBellGate) = 1.0  
op_time(op) = 1.0 # Default case for operations that do not have a defined time

# If true, these ops will be ignored when adding thermal relaxation noise.
is_only_noise(op::PauliNoiseOp) = true
is_only_noise(op::T1NoiseOp) = true
is_only_noise(op::T2NoiseOp) = true
is_only_noise(op) = false # all other operations will be considered non-noise

# "T1, T2 noise: find all 'idle' parts of the QC and fill them with T1 and T2 noise. No noise is applied during a op on that register, but only before and after until the end of the circuit
function noisify_circuit(error::T1T2Noise, circuit; number_registers)
    # tracks how much time has passed on each register since the last operation was applied, so if multiple registers on an incoming operation have the same time, they can receive no T1/T2 noise, but if they have different times, we can add T1/T2 noise to the one that has been idle, effectively simulating a parallel circuit operations with idle noise
    accounted_for_time = zeros(Float64, number_registers) 
    noisy_circuit = [] # to return
    @assert number_registers > 0 # this is assumed
    qubit_list = 1:number_registers # for use in indexing

    # Go through each operation
    for op in circuit 
        # make sure that the operations we care about are real actions on the qubit, not only noise
        # ie, do not add T1T2 noise between other T1T2 Noise! 
        # and ignore any other only-noise operations (ex. PauliNoiseOp)
        if !(is_only_noise(op))
            ## --- Before applying operation ---

            # get time that it would take 
            time_to_elapse = op_time(op)  

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
