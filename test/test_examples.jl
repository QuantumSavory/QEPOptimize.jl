
@testitem "Examples - evaluate" tags=[:examples_plotting] begin
    include("../example/evaluate.jl")
end

@testitem "Examples - optimize" tags=[:examples_plotting] begin
    include("../example/optimize.jl")
end

@testitem "Examples - qasm" begin
    include("../example/QClifford_to_qasm.jl")
end