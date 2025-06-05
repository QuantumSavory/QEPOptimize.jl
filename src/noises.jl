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

"A convenience constructor for unbiased Pauli network noise that results in a Bell pair of given fidelity"
NetworkFidelity(f) = NetworkPauliNoise(f_in_to_pauli(f)...)

function noisify_circuit(n::NetworkPauliNoise, circuit; number_registers)
    initial_noise = [PauliNoiseOp(i, n.px, n.py, n.pz) for i in 1:number_registers]
    return [initial_noise; noisify.((n,), circuit)]
end

noisify(n::NetworkPauliNoise, b::BellMeasure) = NoisyBellMeasureNoisyReset(b, 0, n.px, n.py, n.pz)
noisify(n::NetworkPauliNoise, b::NoisyBellMeasureNoisyReset) = NoisyBellMeasureNoisyReset(b.m, b.p, n.px, n.py, n.pz)


function noisify_circuit(error::MeasurementError, circuit::Vector{Any}; number_registers)
    return map(op -> noisify(error, op), circuit)
end

noisify(error::MeasurementError, b::NoisyBellMeasureNoisyReset) = NoisyBellMeasureNoisyReset(b.m, error.p, b.px, b.py, b.pz)
noisify(error::MeasurementError, b::BellMeasure) = NoisyBellMeasureNoisyReset(b, error.p, 0, 0, 0)

noisify(n::PauliNoise, c::CNOTPerm) = PauliNoiseBellGate(c, n.px, n.py, n.pz) # TODO reusing the QuantumClifford way of making noisy gates would be more consistent here


# TODO more types of noise should be implemented
    # gate_fidelity would turn CNOTPerm gates into gates wrapped into noise
    # T1/T2 noise will add noise that happens even during wait time
    # network_fidelity would turn BellMeasure into NoisyBellMeasureNoisyReset(...)
    # measurement_fidelity would do the same
    # Probably it would be best to define a `noisify(::AbstractNoise, ::AbstractOperation)::AbstractOperation`
    # Or even `noisify_circuit(::AbstractNoise, ::Vector{<:AbstractOperation})` so that we can also cover the T1, T2 noise
