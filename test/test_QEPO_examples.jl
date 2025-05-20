
@testitem "Evaluate example runs" setup=[setup] begin
    include("../example/evaluate.jl")
end

@testitem "Optimize example runs"  begin
    include("../example/optimize.jl")
end

