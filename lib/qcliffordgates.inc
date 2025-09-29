// QuantumClifford gate library for OPENQASM 3.0
// Useful stabilizer transformations in openqasm basis
// OpenQASM 3.0 U gate: exp(i(ϕ+λ)/2) Rz(ϕ)Ry(θ)Rz(λ)
//
// Single-qubit gates present in CNOTPerm:

//  X₁ ⟼ - Y
//  Z₁ ⟼ + X
gate hs q {U(π/2, 0, 1.5*π) q;}

//  X₁ ⟼ + X
//  Z₁ ⟼ - Y
gate plusXminusY q {U(π/2, 1.5*π, π/2) q;}

//  X₁ ⟼ + Z
//  Z₁ ⟼ + Y
gate sh q {U(π/2, π/2, π) q;}

//  X₁ ⟼ + Z
//  Z₁ ⟼ - Y
gate szh q {U(π/2, 1.5*π, π) q;}

//  X₁ ⟼ + X
//  Z₁ ⟼ + Y
gate shs q {U(π/2, π/2, 1.5*π) q;}

//  X₁ ⟼ - Y
//  Z₁ ⟼ + Z
gate sz q {U(0, 0, 1.5*π) q;}

//  X₁ ⟼ + Y
//  Z₁ ⟼ + X
gate hsz q {U(π/2, 0, π/2) q;}

