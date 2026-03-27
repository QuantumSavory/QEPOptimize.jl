// QuantumClifford gate library for OPENQASM 3.0
// Useful stabilizer transformations in openqasm basis
// OpenQASM 3.0 U gate: exp(i(ϕ+λ)/2) Rz(ϕ)Ry(θ)Rz(λ)
//
// Single-qubit gates present in CNOTPerm:

//  X1 -> - Y
//  Z1 -> + X
gate hs q {U(π/2, 0, 1.5*π) q;}

//  X1 -> + X
//  Z1 -> - Y
gate plusXminusY q {U(π/2, 1.5*π, π/2) q;}

//  X1 -> + Z
//  Z1 -> + Y
gate sh q {U(π/2, π/2, π) q;}

//  X1 -> + Z
//  Z1 -> - Y
gate szh q {U(π/2, 1.5*π, π) q;}

//  X1 -> + X
//  Z1 -> + Y
gate shs q {U(π/2, π/2, 1.5*π) q;}

//  X1 -> - Y
//  Z1 -> + Z
gate sz q {U(0, 0, 1.5*π) q;}

//  X1 -> + Y
//  Z1 -> + X
gate hsz q {U(π/2, 0, π/2) q;}

