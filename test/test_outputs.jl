using TestItems

@testitem "qasm bell state preperation" begin
    using QEPOptimize: Φ⁺_qasm, Ψ⁺_qasm, Φ⁻_qasm, Ψ⁻_qasm
    @test Φ⁺_qasm(1,2) == "h q[1];\ncx q[1],q[2];\n"
    @test Φ⁻_qasm(1,2) == "x q[1];\nh q[1];\ncx q[1],q[2];\n"
    @test Ψ⁺_qasm(1,2) == "x q[2];\nh q[1];\ncx q[1],q[2];\n"
    @test Ψ⁻_qasm(1,2) == "x q[1];\nh q[1];\ncx q[1],q[2];\nx q[1];\n"

end

@testitem "to_qasm basic ops" setup=[setup] begin
    using QEPOptimize:to_qasm

    # Test conversion of a single operation CNOTPerm
    qasm_str = to_qasm(setup.test_indiv.ops[1])
    @test occursin("cx q[4], q[8]", qasm_str)

    # Test conversion of the full circuit
    qasm_str_full = to_qasm(setup.test_indiv.ops, setup.number_registers,1)
    @test occursin("cx q[6], q[4]", qasm_str_full)
    # TODO: test measurement when implemented
    # @test occursin("measure q[1] -> c[1]", qasm_str_full)

    # Make sure that header and registers are correct (first four lines)
    @test startswith(qasm_str_full, "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[$(setup.number_registers*2)];creg c[$(setup.number_registers*2)];\n")

end

# TODO
# @testitem "to_qasm CNOTPerm" setup=[setup] begin
#     using QEPOptimize:to_qasm

#     # Test conversion of a single CNOTPerm operation
#     cnot_op = BPGates.CNOTPerm(2, 1, 2, 4)
#     correct_out = "TODO"
#     qasm_str = to_qasm(cnot_op, setup.number_registers)
    
#     # Check that the CNOT operation is equal as the the correct output
#     @test qasm_str == correct_out
# end

# @testitem "to_qasm BellMeasure" setup=[setup] begin
#     using QEPOptimize:to_qasm

#     # Test conversion of a single BellMeasure operation
#     bell_op = BPGates.BellMeasure(1, 1)
#     correct_out = "sdg q[1];\nh q[1];\nmeasure q[1] -> c[1]\n"
#     qasm_str = to_qasm(bell_op)

#     # Check that the BellMeasure operation is equal to the correct output
#     @test qasm_str == correct_out
# end

# TODO t1 t2 errors, noisymeasure turns into BellMeasure, noisy ops turn into clean ops
@testitem "to_qasm ignores noise ops" begin
    using QEPOptimize:to_qasm
    using BPGates:PauliNoise

    # Test conversion of a PauliNoise operation
    noise_op = PauliNoise(0.01/3, 0.01/3, 0.01/3)
    qasm_str = to_qasm(noise_op)

    # Check that the output is an empty string
    @test qasm_str == ""

end

# TODO
@testitem "to_stabilizer" setup=[setup] begin
    using QEPOptimize:to_stabilizer

    output = to_stabilizer(setup.test_ops_small,4)

end
