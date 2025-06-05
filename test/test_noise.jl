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