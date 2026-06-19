# To view the T1T2 errors, Quantikz.displaycircuit requires the BPGates package to export T1NoiseOp and T2NoiseOp, as well as define Quantikz.QuantikzOp for these types. This is done in a fork at https://github.com/xoth42/BPGates.jl

using BPGates

test_ops = [BPGates.CNOTPerm(2, 1, 2, 4), BPGates.CNOTPerm(3, 5, 1, 2), BPGates.CNOTPerm(5, 6, 2, 4), BPGates.CNOTPerm(2, 4, 3, 1), BPGates.CNOTPerm(5, 4, 3, 2), BPGates.CNOTPerm(1, 3, 2, 3), BPGates.CNOTPerm(2, 3, 4, 2), BPGates.CNOTPerm(2, 1, 4, 1), BPGates.BellMeasure(1, 1), BPGates.CNOTPerm(5, 1, 3, 2)]

starting_length = length(test_ops)
registers = 4
pairs = 1:4

using QEPOptimize
using Quantikz


noise_model = QEPOptimize.T1T2Noise(1,1)

# try simplest QC 
test_ops_smallest = [BPGates.CNOTPerm(2, 1, 2, 4),BPGates.BellMeasure(1, 1)]
noisy_circuit_smallest = QEPOptimize.noisify_circuit(noise_model,test_ops_smallest;number_registers=registers)
Quantikz.displaycircuit(test_ops_smallest)
Quantikz.displaycircuit(noisy_circuit_smallest)

# try with small example circuit 
test_ops_small = [BPGates.CNOTPerm(2, 1, 2, 4), BPGates.CNOTPerm(3, 5, 1, 2), BPGates.BellMeasure(1, 1)]
noisy_circuit_small = QEPOptimize.noisify_circuit(noise_model,test_ops_small;number_registers=registers)
Quantikz.displaycircuit(test_ops_small)
Quantikz.displaycircuit(noisy_circuit_small)

# try with larger example circuit
noisy_circuit = QEPOptimize.noisify_circuit(noise_model,test_ops;number_registers=registers)
Quantikz.displaycircuit(test_ops)
Quantikz.displaycircuit(noisy_circuit)

