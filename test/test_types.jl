using TestItems

@testitem "performance averaging" begin
    using QEPOptimize
    using QEPOptimize:Performance

    low = Performance([0,0,1],0,0,0,0,1)
    high = Performance([1,0,0],1,1,1,1,1)
    mid = low + high
    @test mid == Performance([0.5,0,0.5],0.5,0.5,0.5,0.5,2)

    # works with +=
    low = Performance([0,0,1],0,0,0,0,1)
    high = Performance([1,0,0],1,1,1,1,1)

    low += high
    @test low == mid
end
