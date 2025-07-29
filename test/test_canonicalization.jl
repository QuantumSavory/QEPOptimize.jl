using TestItems

@testitem "cleanup_two_measurements"  begin
    using BPGates
    using QEPOptimize:cleanup_two_measurements!

    a = BellMeasure(1,1)
    b = CNOTPerm(1,2,3,4)
    c = BellMeasure(3,0)
    d = BellMeasure(3,0)

    o = [a,b,c,d]

    @test cleanup_two_measurements!(o) == [a,b,c]

    o = [a,b,c,d,d,d,d,d,d,d]
    @test cleanup_two_measurements!(o) == [a,b,c]

    o = [a,b,c,b,d,b,d,b,d,b,d,b,d,b,d,b,d]
    @test cleanup_two_measurements!(o) == [a,b,c,b,d,b,d,b,d,b,d,b,d,b,d,b,d]


    o = [a,d]
    @test cleanup_two_measurements!(o) == [a,d]

    o = [c,c]
    @test cleanup_two_measurements!(o) == [c,c]

end