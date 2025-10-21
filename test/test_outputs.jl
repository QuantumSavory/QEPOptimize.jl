using TestItems

@testmodule QASM begin
    using Conda

    # Add the required Conda channel and install Qiskit if it's not already present.
    Conda.add_channel("conda-forge")
    Conda.add("qiskit")
    Conda.add("qiskit-aer")
    Conda.add("python-dateutil") # deps for qiskit-aer

    # PyCall integration
    using PyCall
    qiskit = pyimport("qiskit")
    AerSimulator = pyimport("qiskit_aer").AerSimulator()

    # helper function to compile/run
    function read(qasm_string) 
        return qiskit.QuantumCircuit.from_qasm_str(qasm_string)
    end
end



@testitem "qasm libs compile" setup=[QASM] begin
    qasm = """
    openqasm 3.0;
    qubit q;
    h q;
    """

    qc = QASM.read(qasm)
    @test isa(qc, PyObject)
    @test Int(qc[:num_qubits]) == 1
end

# @testitem "to_qasm compiles and runs with qiskit_aer - CNOTPerm" begin
#     using QEPOptimize:to_qasm
# 	using BPGates:CNOTPerm,BellMeasure

#     # Try one cnotperm 
#     # qcliff to qasm
#     circ = [CNOTPerm(1,2,3,4)]
#     pairs=4;p_pairs=1;
#     qasm_string = to_qasm(circ,pairs,p_pairs)
#     qasm_string = place_qcliff_lib(qasm_string)


#     # run qasm (run python)
#     # https://quantumcomputing.stackexchange.com/questions/9360/qiskit-in-julia-language
#     using PyCall, Conda

#     Conda.add("pip")
#     Conda.pip("install", "qiskit")
#     qiskit = pyimport("qiskit")

#     # Import Qiskit modules
#     AerSimulator = pyimport("qiskit_aer")

#     # Load circuit from QASM string
#     qc = qiskit.QuantumCircuit.from_qasm_str(qasm_string)

#     # Create Aer simulator
#     simulator = AerSimulator()

#     # Transpile circuit for simulator
#     compiled_circuit = qiskit.transpile(qc, simulator)

#     # Run the circuit with 1000 shots
#     job = simulator.run(compiled_circuit, shots=1000)

#     # Get results
#     result = job.result()
#     counts = result.get_counts(compiled_circuit)

#     println("\nMeasurement Results:")
#     println("Counts: ", counts)
    
# end

@testitem "to_qasm basic ops" setup=[setup] begin
    using QEPOptimize:to_qasm

    # Test conversion of a single operation CNOTPerm
    # This changed: TODO
    # qasm_str = to_qasm(setup.test_indiv.ops[1])
    # @test occursin("cx q[4], q[8]", qasm_str)

    # Test conversion of the full circuit
    qasm_str_full = to_qasm(setup.test_indiv.ops, setup.number_registers,1)
    @test occursin("cx q[6], q[4]", qasm_str_full)
    # TODO: test measurement when implemented
    # @test occursin("measure q[1] -> c[1]", qasm_str_full)

    # Make sure that header and registers are correct (first four lines)
    @test startswith(qasm_str_full, "OPENQASM 3.0;\ninclude \"stdgates.inc\";")

    # Does not fail with some ops/non ops
    Fakeop = 3
    to_qasm(Fakeop)

    # all measure ops contain 'measure'
    using QuantumClifford: sMX,sMZ,sMY
    measure = "measure q"
    @test occursin(measure,to_qasm(sMZ(1)))
    @test occursin(measure,to_qasm(sMX(1)))
    @test occursin(measure,to_qasm(sMY(1)))

end

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
