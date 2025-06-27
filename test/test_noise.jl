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
    using BPGates: affectedqubits
    using BPGates
    # Test with CNOTPerm
    op = BPGates.CNOTPerm(2, 1, 2, 4)
    qubits = affectedqubits(op)
    @test qubits == (2, 4) # correct format

    # Test with BellMeasure
    op = BPGates.BellMeasure(1, 1)
    qubits = affectedqubits(op)
    @test qubits == (1,)

    # Test with operation with no qubits (if such exists)
    struct DummyOp end
    qubits = affectedqubits(DummyOp())
    @test qubits == ()

    # Test with NoisyBellMeasureNoisyReset
    op = BPGates.NoisyBellMeasureNoisyReset(BPGates.BellMeasure(1, 1), 0.1, 0.01, 0.01, 0.01)
    qubits = affectedqubits(op)
    @test qubits == (1,)

    # Test with PauliNoiseBellGate
    op = BPGates.PauliNoiseBellGate(BPGates.CNOTPerm(2, 1, 2, 4), 0.01, 0.01, 0.01)
    qubits = affectedqubits(op)
    @test qubits == (2, 4)
    
end

# @testitem "T1T2 Noise obeys noise equality" begin
#     using QEPOptimize: thermal_relaxation_error_rate
#     t1 = 1;t2 = 1;gate_time = 1
#     λ₁,λ₂ = thermal_relaxation_error_rate(t1,t2,gate_time)
#     doubleλ₁,doubleλ₂ = thermal_relaxation_error_rate(t1,t2,gate_time * 2)

#     # Must obey lambda_2x = 2*lambda_x*(1-lambda_x)
#     # λ₁
#     @test doubleλ₁ == 2 * λ₁ * (1 - λ₁)

#     # λ₂
#     @test doubleλ₂ == 2 * λ₂ * (1 - λ₂)

# end

@testitem "Add T1T2 noise" setup=[noise_setup] begin
    using QEPOptimize: T1T2Noise, noisify_circuit
    using BPGates
    noise_model = T1T2Noise(1,1)

    noisy_circuit = noisify_circuit(noise_model,noise_setup.test_ops;number_registers=noise_setup.registers)

    # Check that the length of the noisy circuit is greater than the original
    @test length(noisy_circuit) > length(noise_setup.test_ops)

    # Check that the noisy circuit contains the correct amount of T1T2Noise operations
    t1_ops = filter(op -> isa(op, BPGates.T1NoiseOp), noisy_circuit)
    @test length(t1_ops) == 9
    t2_ops = filter(op -> isa(op, BPGates.T2NoiseOp), noisy_circuit)
    @test length(t2_ops) == 9

    # test with small circuit example
    noisy_circuit_small = noisify_circuit(noise_model,noise_setup.test_ops_small;number_registers=noise_setup.registers)
    # Check that the length of the noisy circuit is greater than the original
    @test length(noisy_circuit_small) > length(noise_setup.test_ops_small)
    # Check that the noisy circuit contains the correct amount of T1T2Noise operations
    t1_ops_small = filter(op -> isa(op, BPGates.T1NoiseOp), noisy_circuit_small)
    @test length(t1_ops_small) == 1
    t2_ops_small = filter(op -> isa(op, BPGates.T2NoiseOp), noisy_circuit_small)
    @test length(t2_ops_small) == 1

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
    
    # Check default times for different operations
    @test op_time(default_times.op_times, BPGates.CNOTPerm(1, 2, 3, 4)) == 1.0
    @test op_time(default_times.op_times, BPGates.BellMeasure(1, 2)) == 1.0
    @test op_time(default_times.op_times, BPGates.NoisyBellMeasureNoisyReset(BPGates.BellMeasure(1, 2), 0.1, 0.01, 0.01, 0.01)) == 1.0
    @test op_time(default_times.op_times, BPGates.PauliNoiseBellGate(BPGates.CNOTPerm(1, 2, 3, 4), 0.1, 0.01, 0.01)) == 1.0
    
    # Noise operations should have zero time by default
    @test op_time(default_times.op_times, BPGates.T1NoiseOp(1, 0.1)) == 0.0
    @test op_time(default_times.op_times, BPGates.T2NoiseOp(1, 0.1)) == 0.0
    @test op_time(default_times.op_times, BPGates.PauliNoiseOp(1, 0.1, 0.1, 0.1)) == 0.0

end
    
