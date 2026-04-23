using QEPOptimize:to_qasm
using BPGates:CNOTPerm,BellMeasure


my_circuit = [
    CNOTPerm(2, 1, 2, 4),
    CNOTPerm(3, 5, 1, 2),
    CNOTPerm(5, 6, 2, 4),
    CNOTPerm(2, 4, 3, 1), 
    CNOTPerm(5, 4, 3, 2), 
    CNOTPerm(1, 3, 2, 3),
    CNOTPerm(2, 3, 4, 2), 
    CNOTPerm(2, 1, 4, 1), 
    BellMeasure(1, 1), 
    CNOTPerm(5, 1, 3, 2)
]

pairs = 4
purified_pairs = 1

output = to_qasm(my_circuit,pairs,purified_pairs)
print(output)
