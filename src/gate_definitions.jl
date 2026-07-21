using QuantumClifford

# Common gates translation (SparseGate)
# Note: in this code we use 'p/P' as shorthand for the phase gate or sqrt Z. In qasm, this gate is called 's'. We will try to use s when possible to reduce confusion. 

# Custom single-qubit gate definitions
const tHS = C"-Y 
               X"
# X₁ ⟼ - Y
# Z₁ ⟼ + X
# [	0.707	0.707im
# 	0.707	-0.707im	]
# U(0.500π, 0, 1.500π)

const tplusXminusY = C"X
                     -Y"
# X₁ ⟼ + X
# Z₁ ⟼ - Y
# [	0.707	-0.707im
# 	-0.707im	0.707	]
# U(0.500π, 1.500π, 0.500π)   

const tSH = C"Z
              Y"
# X₁ ⟼ + Z
# Z₁ ⟼ + Y
# [	0.707	0.707
# 	0.707im	-0.707im	]
# U(0.500π, 0.500π, 1.000π)

const tHSZ = C"Y
               X"         
# X₁ ⟼ + Y
# Z₁ ⟼ + X
# [	0.707	-0.707im
# 	0.707	0.707im	]
# U(0.500π, 0, 0.500π)

const tSZH = C"Z
                -Y"
# X₁ ⟼ + Z
# Z₁ ⟼ - Y
# [	0.707	0.707
# 	-0.707im	0.707im	]
# U(0.500π, 1.500π, 1.000π)

const tSZ = C"-Y    
                Z"
# X₁ ⟼ - Y
# Z₁ ⟼ + Z
# [	1.000	0
# 	0	-1.000im	]
# U(0, 0, 1.500π)

const tSHS =  C"X 
                Y"
# X₁ ⟼ + X
# Z₁ ⟼ + Y
# [	0.707	0.707im
# 	0.707im	0.707	]
# U(0.500π, 0.500π, 1.500π)

# Gate mapping dictionaries
# common single-qubit sparsegates to qasm, only gate name
const single_SparseGate_to_qasm = Dict(
    tHS =>"hs",
    tplusXminusY =>"plusXminusY",
    tSH =>"sh",
    tHSZ =>"hsz",
    tSZH =>"szh",
    tSZ =>"sz",
    tSHS =>"shs",
    tHadamard => "h",
    tPhase => "s"
)

const double_SparseGate_to_qasm = Dict(
    tCNOT => "cx",
    tSWAP => "swap",
    tCPHASE => "cz"
)
