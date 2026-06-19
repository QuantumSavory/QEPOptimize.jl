using TestItems

@testmodule noise_setup begin
    using BPGates

    test_ops = [BPGates.CNOTPerm(2, 1, 2, 4), BPGates.CNOTPerm(3, 5, 1, 2), BPGates.CNOTPerm(5, 6, 2, 4), BPGates.CNOTPerm(2, 4, 3, 1), BPGates.CNOTPerm(5, 4, 3, 2), BPGates.CNOTPerm(1, 3, 2, 3), BPGates.CNOTPerm(2, 3, 4, 2), BPGates.CNOTPerm(2, 1, 4, 1), BPGates.BellMeasure(1, 1), BPGates.CNOTPerm(5, 1, 3, 2)]

    test_ops_small = [BPGates.CNOTPerm(2, 1, 2, 4),BPGates.BellMeasure(1, 1)]

    starting_length = length(test_ops)
    registers = 4
    pairs = 1:4

end

@testitem "noisify - Measurement Error" setup=[noise_setup] begin
    using QEPOptimize: noisify_circuit, MeasurementError
    using BPGates
    # Test with zero error rate (should not add any operations)
    noise_model = MeasurementError(0.0)
    noisy_ops = noisify_circuit(noise_model,noise_setup.test_ops;number_registers=noise_setup.registers)
    @test length(filter(op -> isa(op,BPGates.NoisyBellMeasureNoisyReset) ? op.p != 0 : false,noisy_ops)) == 0
    
    # Test with non-zero error rate
    noise_model = MeasurementError(0.1)
    noisy_ops = noisify_circuit(noise_model,noise_setup.test_ops;number_registers=noise_setup.registers)
    @test length(noisy_ops) == noise_setup.starting_length
    
    # Test that only measurement operations are affected
    # Find all BellMeasure operations in original circuit
    measure_indices = findall(op -> isa(op, BPGates.BellMeasure), noise_setup.test_ops)
    
    # Check that noise was applied around these operations
    for idx in measure_indices
        @test isa(noisy_ops[idx], BPGates.NoisyBellMeasureNoisyReset)
        @test noisy_ops[idx].p == noise_model.p
    end
    
    # Test with 100% error rate
    noise_model = MeasurementError(1.0)
    noisy_ops = noisify_circuit(noise_model,noise_setup.test_ops;number_registers=noise_setup.registers)
    @test length(noisy_ops) == noise_setup.starting_length

    # simulation should fail completely
    simulations = 100
    for _ in 1:simulations
            
        purified_state, res = BPGates.mctrajectory!(BPGates.BellState(noise_setup.registers), noisy_ops)
            # If the circuit execution was reported as 'successful'
        @test res != BPGates.continue_stat
    end
    
end

@testitem "affectedqubits" begin
    using Quantikz: affectedqubits
    using BPGates
    # Test with CNOTPerm
    op = BPGates.CNOTPerm(2, 1, 2, 4)
    qubits = affectedqubits(op)
    @test qubits == [2, 4]

    # Test with BellMeasure
    op = BPGates.BellMeasure(1, 1)
    qubits = affectedqubits(op)
    @test qubits == [1]

    # Test with NoisyBellMeasureNoisyReset
    op = BPGates.NoisyBellMeasureNoisyReset(BPGates.BellMeasure(1, 1), 0.1, 0.01, 0.01, 0.01)
    qubits = affectedqubits(op)
    @test qubits == [1]

    # Test with PauliNoiseBellGate
    op = BPGates.PauliNoiseBellGate(BPGates.CNOTPerm(2, 1, 2, 4), 0.01, 0.01, 0.01)
    qubits = affectedqubits(op)
    @test qubits == [2, 4]
    
end

@testitem "T1T2 Noise obeys noise equality" begin
    using QEPOptimize: thermal_relaxation_error_rate
    t1 = 1;t2 = 1;gate_time = 1
    λ₁,λ₂ = thermal_relaxation_error_rate(t1,t2,gate_time)
    doubleλ₁,doubleλ₂ = thermal_relaxation_error_rate(t1,t2,gate_time * 2)

    # evolving for 2t is equivalent to two steps of t, so S(2t)=S(t)^2 ⇒ λ(2t)=1-(1-λ(t))^2
    # λ₁
    @test doubleλ₁ == 1 - (1 - λ₁)^2

    # λ₂
    @test doubleλ₂ == 1 - (1 - λ₂)^2

end

"""
T1 noise inaccurate with low t1, only use >100 for roughly correct t1 noise. T2 is fine. There is some nuance here. The T1 Noise model uses a markov chain made from the populations of analytic results of T1 decay with amplitude damping kraus operators:
 `|0⟩⟨0| + √(1-λ) |1⟩⟨1|` and `√λ |0⟩⟨1|`
It acts on bell states, and deriving the full results we can look at one example:
Input state: 00 or |Φ+⟩, Amplitude damped Output density:
┌                                                       ┐
│ 0.5λ₁² - λ₁ + 1  0.5λ₁       0            0           │
│ 0.5λ₁            0.5λ₁²      0            0           │
│ 0                0           0.5λ₁(1-λ₁)  0           │
│ 0                0           0            0.5λ₁(1-λ₁) │
└                                                       ┘
There are non-zero coherences (0.5λ₁) which are not modeled in the chain (only uses populations). These are damped out normally with Tᵩ (See Manenti & Motta, 8.4-6). If T1 is very small, these coherences are large and the model is inaccurate. If T1 is large, these coherences are small and the model is accurate enough.

The T2 Noise model (actually Tᵩ - isolated phase damping) works perfectly, because there are no coherences. Example:
Input state: 00 or |Φ+⟩, Phase damped output density:
┌                                                   ┐
│ 0.5λ₂² - λ₂ + 1    0                0         0   │
│ 0                  λ₂(1 - 0.5λ₂)    0         0   │
│ 0                  0                0         0   │
│ 0                  0                0         0   │
└                                                   ┘
"""
@testitem "two 1-sec T1 noise has same results as one 2-sec T1 noise" begin
    using QEPOptimize: thermal_relaxation_error_rate,Individual,calculate_performance!
    using BPGates: T1NoiseOp

    two_noise_ops_T1 = Individual()
    one_noise_op_T1 = Individual()

    # T1 must fail on bad error cases, but it will work on normal T1 values
    t1 = 100 
    t2 = 150 
    gate_time = 0.2 # very long wait time in μs 
    λ₁, λ₂ = thermal_relaxation_error_rate(t1, t2, gate_time)

    # Create two T1 noise ops
    two_noise_ops_T1.ops = [T1NoiseOp(1, λ₁), T1NoiseOp(1, λ₁)]
    # Create one T1 noise op with double time
    double_λ₁, _ = thermal_relaxation_error_rate(t1, t2, gate_time * 2)
    one_noise_op_T1.ops = [T1NoiseOp(1, double_λ₁)]

    sims = 100000
    one_sec_T1_performance = calculate_performance!(two_noise_ops_T1; num_simulations=sims, noises=[])
    two_sec_T1_performance = calculate_performance!(one_noise_op_T1; num_simulations=sims, noises=[])

    @test isapprox(one_sec_T1_performance.logical_qubit_fidelity, two_sec_T1_performance.logical_qubit_fidelity; atol=5/sqrt(sims))

end

@testitem "two 1-sec T2 noise has same results as one 2-sec T2 noise" begin
    using QEPOptimize: thermal_relaxation_error_rate,Individual,calculate_performance!
    using BPGates: T2NoiseOp

    t1 = 1
    t2 = 1
    gate_time = 1
    one_sec_λ₁, one_sec_λ₂ = thermal_relaxation_error_rate(t1, t2, gate_time)
    two_sec_λ₁, two_sec_λ₂ = thermal_relaxation_error_rate(t1, t2, gate_time*2)
    
    two_noise_ops_T2 = Individual()
    two_noise_ops_T2.ops = [T2NoiseOp(1,one_sec_λ₂),T2NoiseOp(1,one_sec_λ₂)]

    one_noise_op_T2 = Individual()
    one_noise_op_T2.ops = [T2NoiseOp(1,two_sec_λ₂)]

    sims = 100000
    one_sec_T2_performance = calculate_performance!(two_noise_ops_T2; num_simulations=sims,noises=[])
    two_sec_T2_performance = calculate_performance!(one_noise_op_T2; num_simulations=sims,noises=[])
    
    tolerance = 5/sqrt(sims)
    @test isapprox(one_sec_T2_performance.logical_qubit_fidelity, two_sec_T2_performance.logical_qubit_fidelity;atol=tolerance)
end


@testitem "Add T1T2 noise" setup=[noise_setup] begin
    using QEPOptimize: T1T2Noise, noisify_circuit
    using BPGates
    noise_model = T1T2Noise(1,1)

    test_ops = noise_setup.test_ops
    registers = noise_setup.registers
    noisy_circuit = noisify_circuit(noise_model,test_ops;number_registers=registers)

    # Check that the length of the noisy circuit is greater than the original
    @test length(noisy_circuit) > length(test_ops)

    # Check that the noisy circuit contains the correct amount of T1T2Noise operations
    t1_ops = filter(op -> isa(op, BPGates.T1NoiseOp), noisy_circuit)
    @test length(t1_ops) == 9
    t2_ops = filter(op -> isa(op, BPGates.T2NoiseOp), noisy_circuit)
    @test length(t2_ops) == 9

    # test with small circuit example
    test_ops_small = noise_setup.test_ops_small
    noisy_circuit_small = noisify_circuit(noise_model,test_ops_small;number_registers=registers)
    # Check that the length of the noisy circuit is greater than the original
    @test length(noisy_circuit_small) > length(test_ops_small)
    # Check that the noisy circuit contains the correct amount of T1T2Noise operations
    t1_ops_small = filter(op -> isa(op, BPGates.T1NoiseOp), noisy_circuit_small)
    @test length(t1_ops_small) == 3
    t2_ops_small = filter(op -> isa(op, BPGates.T2NoiseOp), noisy_circuit_small)
    @test length(t2_ops_small) == 3
    print(t1_ops_small)
    print(noisy_circuit_small)

end

@testitem "op_times for T1T2 Noise" begin
    using QEPOptimize: T1T2Noise,op_time
    using BPGates
    
    # Test custom op_times
    times = Dict(
        BellMeasure => 1.0,
        CNOTPerm => 2.0,
        PauliNoiseBellGate => 3.0,
        NoisyBellMeasureNoisyReset => 4.0
    )
    
    # Noise ops: T1NoiseOp, T2NoiseOp, PauliNoiseOp
    # Non-noise ops: CNOTPerm, BellMeasure, NoisyBellMeasureNoisyReset, PauliNoiseBellGate
    # Noise:
    @test op_time(times, BPGates.T1NoiseOp(1, 0.1)) == 0.0
    @test op_time(times, BPGates.T2NoiseOp(1, 0.1)) == 0.0
    @test op_time(times, BPGates.PauliNoiseOp(1, 0.1, 0.1, 0.1)) == 0.0
    @test op_time(times, BPGates.CNOTPerm(1, 2, 3, 4)) == 2.0
    @test op_time(times, BPGates.BellMeasure(1, 2)) == 1.0
    @test op_time(times, BPGates.NoisyBellMeasureNoisyReset(BPGates.BellMeasure(1, 2), 0.1, 0.01, 0.01, 0.01)) == 4.0
    @test op_time(times, BPGates.PauliNoiseBellGate(BPGates.CNOTPerm(1, 2, 3, 4), 0.1, 0.01, 0.01)) == 3.0


    # Test default op_times
    default_times = T1T2Noise(1.0, 2.0)
    
    # Noise operations should have zero time by default
    @test op_time(default_times.op_times, BPGates.T1NoiseOp(1, 0.1)) == 0.0
    @test op_time(default_times.op_times, BPGates.T2NoiseOp(1, 0.1)) == 0.0
    @test op_time(default_times.op_times, BPGates.PauliNoiseOp(1, 0.1, 0.1, 0.1)) == 0.0

end
    
