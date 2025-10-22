using TestItems

@testmodule QASM begin
    using Conda

    @doc """
    Test utilities for QASM output validation and integration testing with Qiskit.
    Provides setup functions and helpers for compiling/running QASM code.
    """
    function setup_qiskit_conda()
        # Add the required Conda channel and install Qiskit if it's not already present.
        if !("conda-forge" in Conda.channels())
            Conda.add_channel("conda-forge")
        end
        
        # Check if packages are already installed before adding them
        installed_pkgs = Conda.parseconda(`list`)
        installed_pkg_names = [pkg["name"] for pkg in installed_pkgs]
        pkgs = ["qiskit","qiskit-aer","python-dateutil","qiskit-qasm3-import"]
        for name in pkgs
            if !(name in installed_pkg_names)
                Conda.add(name)
            end
        end
    end
    
    setup_qiskit_conda()
    # PyCall integration
    using PyCall
    qiskit = pyimport("qiskit")
    sim = pyimport("qiskit_aer").AerSimulator()
    qasm3_parse = pyimport("qiskit_qasm3_import").parse


    # helper function to compile/run
    function compile_qasm(qasm_string) 
        # Handle OPENQASM 3.0 vs openqasm 3.0 case sensitivity
        qasm_string_caps = replace(qasm_string, r"(?i)^openqasm" => "OPENQASM")
        # For qasm 2.0
        # return qiskit.QuantumCircuit.from_qasm_str(QASM.expand_includes(qasm_string))
        # QASM 3.0 needs this wrapper 
        return qasm3_parse(qasm_string_caps)
    end
    
    # helper function to replace includes with actual file contents
    function expand_includes(qasm_string)
        # Get the project root directory (two levels up from test directory)
        project_root = dirname(@__DIR__)
        qasm_string = replace(qasm_string, "include \"stdgates.inc\";" => read(joinpath(project_root, "lib", "stdgates.inc"), String))
        qasm_string = replace(qasm_string, "include \"qcliffordgates.inc\";" => read(joinpath(project_root, "lib", "qcliffordgates.inc"), String))
        return qasm_string
    end
end

@testitem "qasm libs compile" setup=[QASM] tags=[:qasm, :integration] begin
    using PyCall
    qasm = """
    openqasm 3.0;
    include \"stdgates.inc\";
    qubit q;
    h q;
    """

    qc = QASM.compile_qasm(qasm)
    @test isa(qc, PyObject)
    @test Int(qc[:num_qubits]) == 1
end

@testitem "CNOTPerm to_qasm compiles and runs with qiskit_aer" setup=[QASM] tags=[:qasm, :integration] begin
    using QEPOptimize:to_qasm
	using BPGates:CNOTPerm,BellMeasure, sMZ

    # Try one cnotperm 
    # qcliff to qasm
    circ = [CNOTPerm(1,2,3,4),sMZ(1),sMZ(2),sMZ(3),sMZ(4)]
    pairs=4;p_pairs=1;
    qasm_string = to_qasm(circ,pairs,p_pairs)
    qasm_string = QASM.expand_includes(qasm_string)


    # run qasm (run python)
    # https://quantumcomputing.stackexchange.com/questions/9360/qiskit-in-julia-language
    @info qasm_string
    # Load circuit from QASM string
    qc = QASM.compile_qasm(qasm_string)

    compiled_circuit =  QASM.qiskit.transpile(qc, QASM.sim)

    # Run the circuit with 1000 shots
    job = QASM.sim.run(compiled_circuit, shots=1000)

    # Get results
    result = job.result()
    counts = result.get_counts(compiled_circuit)

    println("\nMeasurement Results:")
    println("Counts: ", counts)
end
    
@testitem "stdgates.inc compiles and runs with qiskit_aer" setup=[QASM] tags=[:qasm, :integration] begin
    using QEPOptimize: to_qasm
    using QuantumClifford: sCNOT
    using BPGates: sMZ

    # Test basic stdgates.inc operations
    circ = [sCNOT(1, 2),sMZ(1),sMZ(2)]
    pairs = 2
    p_pairs = 1
    qasm_string = to_qasm(circ, pairs, p_pairs)
    qasm_string = QASM.expand_includes(qasm_string)

    # Load circuit from QASM string
    @info qasm_string
    qc = QASM.compile_qasm(qasm_string)
    compiled_circuit = QASM.qiskit.transpile(qc, QASM.sim)

    # Run the circuit with 1000 shots
    job = QASM.sim.run(compiled_circuit, shots=1000)

    # Get results
    result = job.result()
    counts = result.get_counts(compiled_circuit)

    println("\nstdgates.inc CNOT Results:")
    println("Counts: ", counts)

    # Verify the circuit was created and has expected properties
    @test Int(qc[:num_qubits]) >= 2
    @test haskey(counts, "0000") || haskey(counts, "0011") # Bell state outcomes
end

@testitem "qcliffordgates.inc compiles and runs with qiskit_aer" setup=[QASM] tags=[:qasm, :integration] begin
    using QEPOptimize: to_qasm
    using QuantumClifford: SparseGate
    using QEPOptimize: tHS,sMZ  # Import custom gate from gate_definitions.jl

    # Test custom qcliffordgates.inc operations  
    circ = [SparseGate(tHS, (1,)),sMZ(1)]  # HS gate on qubit 1
    pairs = 2
    p_pairs = 1
    qasm_string = to_qasm(circ, pairs, p_pairs)
    qasm_string = QASM.expand_includes(qasm_string)

    # Load circuit from QASM string
    qc = QASM.compile_qasm(qasm_string)
    compiled_circuit = QASM.qiskit.transpile(qc, QASM.sim)

    # Run the circuit with 1000 shots
    job = QASM.sim.run(compiled_circuit, shots=1000)

    # Get results
    result = job.result()
    counts = result.get_counts(compiled_circuit)

    println("\nqcliffordgates.inc HS gate Results:")
    println("Counts: ", counts)

    # Verify the circuit was created
    @test Int(qc[:num_qubits]) >= 2
end

@testitem "BellMeasure compiles and runs with qiskit_aer" setup=[QASM] tags=[:qasm, :integration] begin
    using QEPOptimize: to_qasm
    using BPGates: BellMeasure

    # Test BellMeasure operation (measurement in X basis)
    circ = [BellMeasure(1, 1)]  # Bell measure on pair 1, X measurement
    pairs = 2
    p_pairs = 1
    qasm_string = to_qasm(circ, pairs, p_pairs)
    qasm_string = QASM.expand_includes(qasm_string)

    # Load circuit from QASM string
    qc = QASM.compile_qasm(qasm_string)
    compiled_circuit = QASM.qiskit.transpile(qc, QASM.sim)

    # Run the circuit with 1000 shots
    job = QASM.sim.run(compiled_circuit, shots=1000)

    # Get results
    result = job.result()
    counts = result.get_counts(compiled_circuit)

    println("\nBellMeasure Results:")
    println("Counts: ", counts)

    # Verify the circuit has measurements
    @test Int(qc[:num_qubits]) >= 2
    @test occursin("measure", qasm_string)  # Should contain measurement operations
end

@testitem "Quantum measurements compiles and runs with qiskit_aer" setup=[QASM] tags=[:qasm, :integration] begin
    using QEPOptimize: to_qasm
    using QuantumClifford: sMX, sMY, sMZ

    # Test different measurement bases
    circ = [sMX(1), sMY(2), sMZ(3)]  # X, Y, Z measurements on different qubits
    pairs = 3
    p_pairs = 1
    qasm_string = to_qasm(circ, pairs, p_pairs)
    qasm_string = QASM.expand_includes(qasm_string)

    # Load circuit from QASM string
    qc = QASM.compile_qasm(qasm_string)
    compiled_circuit = QASM.qiskit.transpile(qc, QASM.sim)

    # Run the circuit with 1000 shots
    job = QASM.sim.run(compiled_circuit, shots=1000)

    # Get results
    result = job.result()
    counts = result.get_counts(compiled_circuit)

    println("\nQuantum Measurements (X,Y,Z) Results:")
    println("Counts: ", counts)

    # Verify measurements are present
    @test occursin("measure", qasm_string)
    @test occursin("h q[0]", qasm_string)  # X measurement requires hadamard
    @test occursin("sdg q[2]", qasm_string)  # Y measurement requires S†
    @test Int(qc[:num_qubits]) >= 3
end

@testitem "Mixed circuit compiles and runs with qiskit_aer" setup=[QASM] tags=[:qasm, :integration] begin
    using QEPOptimize: to_qasm
    using QuantumClifford: sCNOT, SparseGate, sMZ
    using QEPOptimize: tSH  # Import custom gate

    # Test mixed circuit with different gate types
    circ = [
        SparseGate(tSH, (1,)),  # Custom gate on qubit 1
        sCNOT(2, 3),           # CNOT between qubits 2 and 3
        sMZ(1)                 # Z measurement on qubit 1
    ]
    pairs = 3
    p_pairs = 2
    qasm_string = to_qasm(circ, pairs, p_pairs)
    qasm_string = QASM.expand_includes(qasm_string)

    # Load circuit from QASM string
    @info qasm_string
    qc = QASM.compile_qasm(qasm_string)
    compiled_circuit = QASM.qiskit.transpile(qc, QASM.sim)

    # Run the circuit with 1000 shots
    job = QASM.sim.run(compiled_circuit, shots=1000)

    # Get results
    result = job.result()
    counts = result.get_counts(compiled_circuit)

    println("\nMixed Circuit Results:")
    println("Counts: ", counts)

    # Verify circuit properties
    @test occursin("sh q[0]", qasm_string)     # Custom SH gate
    @test occursin("cx q[2], q[4]", qasm_string) # CNOT gate
    @test occursin("measure q[0]", qasm_string) # Measurement
    @test Int(qc[:num_qubits]) >= 3
end


@testitem "to_qasm basic ops" setup=[setup] tags=[:qasm] begin
    using QEPOptimize:to_qasm
    using BPGates: CNOTPerm, BellMeasure
    # Test conversion of a single operation CNOTPerm
    # This changed: TODO
    # qasm_str = to_qasm(setup.test_indiv.ops[1])
    # @test occursin("cx q[4], q[8]", qasm_str)

    # Test conversion of the full circuit
    qasm_str_full = to_qasm([CNOTPerm(1,2,3,4),BellMeasure(1,4)], setup.number_registers,1)
    @test occursin("cx q[4], q[6]", qasm_str_full)
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
@testitem "to_qasm ignores noise ops" tags=[:qasm] begin
    using QEPOptimize:to_qasm
    using BPGates:PauliNoise

    # Test conversion of a PauliNoise operation
    noise_op = PauliNoise(0.01/3, 0.01/3, 0.01/3)
    qasm_str = to_qasm(noise_op)

    # Check that the output is an empty string
    @test qasm_str == ""

end

@testitem "stdgates.inc library" setup=[QASM] tags=[:qasm, :integration] begin
    using PyCall
    # Test basic gates from stdgates.inc
    qasm = """
    OPENQASM 3.0;
    include "stdgates.inc";
    qubit[2] q;
    h q[0];
    x q[1];
    cx q[0], q[1];
    """
    
    qc = QASM.compile_qasm(qasm)
    @test isa(qc, PyObject)
    @test Int(qc[:num_qubits]) == 2
end

@testitem "qcliffordgates.inc library compiles" setup=[QASM] tags=[:qasm, :integration] begin
    using PyCall
    # Test custom gates from qcliffordgates.inc
    qasm = """
    OPENQASM 3.0;
    include "stdgates.inc";
    include "qcliffordgates.inc";
    qubit q;
    hs q;
    plusXminusY q;
    sh q;
    """
    qasm = QASM.expand_includes(qasm)
    qc = QASM.compile_qasm(qasm)
    @test isa(qc, PyObject)
    @test Int(qc[:num_qubits]) == 1
end

@testitem "stdgates.inc single qubit gates" setup=[QASM] tags=[:qasm, :integration] begin
    using PyCall
    # Test individual single-qubit gates from stdgates.inc
    gates = ["h", "x", "y", "z", "s", "sdg", "t", "tdg", "sx"]
    
    for gate in gates
        qasm = """
        OPENQASM 3.0;
        include "stdgates.inc";
        qubit q;
        $gate q;
        """
        
        qc = QASM.compile_qasm(qasm)
        @test isa(qc, PyObject)
        @test Int(qc[:num_qubits]) == 1
    end
end

@testitem "stdgates.inc two qubit gates" setup=[QASM] tags=[:qasm, :integration] begin
    using PyCall
    # Test two-qubit gates from stdgates.inc
    gates = ["cx", "cy", "cz", "swap"]
    
    for gate in gates
        qasm = """
        OPENQASM 3.0;
        include "stdgates.inc";
        qubit[2] q;
        $gate q[0], q[1];
        """
        
        qc = QASM.compile_qasm(qasm)
        @test isa(qc, PyObject)
        @test Int(qc[:num_qubits]) == 2
    end
end

@testitem "qcliffordgates.inc custom gates" setup=[QASM] tags=[:qasm, :integration] begin
    using PyCall
    # Test each custom gate from qcliffordgates.inc
    gates = ["hs", "plusXminusY", "sh", "szh", "shs", "sz", "hsz"]
    qasm_start = """
    OPENQASM 3.0;
    include "stdgates.inc";
    include "qcliffordgates.inc";
    qubit q;
    """
    qasm_start =  QASM.expand_includes(qasm_start)

    for gate in gates
        qasm = qasm_start * "$gate q;"
        
        
        qc = QASM.compile_qasm(qasm)
        @test isa(qc, PyObject)
        @test Int(qc[:num_qubits]) == 1
    end
end

@testitem "expand_includes function" setup=[QASM] tags=[:qasm, :integration] begin
    # Test the expand_includes function directly
    qasm_with_includes = """
    OPENQASM 3.0;
    include "stdgates.inc";
    include "qcliffordgates.inc";
    qubit q;
    h q;
    """
    
    expanded_qasm = QASM.expand_includes(qasm_with_includes)
    
    # Check that includes are replaced with actual content
    @test !occursin("include \"stdgates.inc\";", expanded_qasm)
    @test !occursin("include \"qcliffordgates.inc\";", expanded_qasm)
    
    # Check that actual gate definitions are present
    @test occursin("gate h a", expanded_qasm)
    @test occursin("gate hs q", expanded_qasm)
    
    # Check that the original circuit content is preserved
    @test occursin("qubit q;", expanded_qasm)
    @test occursin("h q;", expanded_qasm)
end

@testitem "replace_special_chars function" tags=[:qasm] begin
    using QEPOptimize: replace_special_chars
    
    # Test with superscripts
    test_str_super = "h q[0]; // Apply H gate to qubit q₁"
    result_super = replace_special_chars(test_str_super)
    @test result_super == "h q[0]; // Apply H gate to qubit q1"
    
    # Test with subscripts
    test_str_sub = "cx q[0], q[1]; // CNOT from q¹ to q²"
    result_sub = replace_special_chars(test_str_sub)
    @test result_sub == "cx q[0], q[1]; // CNOT from q1 to q2"
    
    # Test with both
    test_str_both = "measure q₁ -> c¹;"
    result_both = replace_special_chars(test_str_both)
    @test result_both == "measure q1 -> c1;"
    
    # Test with no special characters
    test_str_normal = "h q[0]; cx q[0], q[1];"
    result_normal = replace_special_chars(test_str_normal)
    @test result_normal == test_str_normal
end
