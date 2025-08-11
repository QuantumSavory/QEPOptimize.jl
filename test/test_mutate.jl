using TestItems

@testmodule setup begin
    using BPGates
    using QEPOptimize:Individual

    test_ops = [BPGates.CNOTPerm(2, 1, 2, 4), BPGates.CNOTPerm(3, 5, 1, 2), BPGates.CNOTPerm(5, 6, 2, 4), BPGates.CNOTPerm(2, 4, 3, 1), BPGates.CNOTPerm(5, 4, 3, 2), BPGates.CNOTPerm(1, 3, 2, 3), BPGates.CNOTPerm(2, 3, 4, 2), BPGates.CNOTPerm(2, 1, 4, 1), BPGates.BellMeasure(1, 1), BPGates.CNOTPerm(5, 1, 3, 2)]

    test_indiv = Individual(test_ops)

    test_ops_small = [BPGates.CNOTPerm(2, 1, 2, 4),BPGates.BellMeasure(1, 1)]
    test_indiv_small = Individual(test_ops_small)
    
    number_registers = 4

    starting_length = length(test_ops)
    empty = Individual()

    pairs = 1:4

end

@testitem "drop_op" setup=[setup] begin
    using QEPOptimize:drop_op

    # normal drop
    @test length(drop_op(setup.test_indiv).ops) == (setup.starting_length - 1) 

    # handle empty ops
    @test setup.empty.ops == drop_op(setup.empty).ops

end

@testitem "gain_op" setup=[setup] begin
    using QEPOptimize:gain_op

    # normal gain
    @test length(gain_op(setup.test_indiv;valid_pairs=setup.pairs).ops) == (setup.starting_length + 1)

    # test with empty individual
    @test length(gain_op(setup.empty;valid_pairs=setup.pairs).ops) == 1
end


@testitem "rand_op" setup=[setup] begin
    using QEPOptimize:rand_op
    using BPGates
    
    # Check that rand_op returns a valid operation
    for _ in 1:100 
        ops = rand_op(setup.pairs)
        @test (ops isa BPGates.CNOTPerm ) | (ops isa BPGates.BellOp ) | (ops isa BPGates.BellMeasure)
    end
    # Test multiple calls to ensure variety
    ops = [rand_op(setup.pairs) for _ in 1:10]
    
    # Test adding the random ops to an individual
    ops_indiv = Individual(ops)
    @test length(ops_indiv.ops) == 10

    # At least some operations should be different (very low probability of all being the same)
    @test length(unique(typeof.(ops))) > 1

end

@testitem "mutate" setup=[setup] begin
    using QEPOptimize:mutate
    
    # Test mutation produces a different individual
    mutated = mutate(setup.test_indiv)
    @test mutated isa Individual
    @test mutated.ops != setup.test_indiv.ops
    
    # Test mutation on empty individual
    mutated_empty = mutate(setup.empty)
    @test mutated_empty isa Individual
    # Empty individual should not gain an operation when mutated
    @test isempty(mutated_empty.ops)
    
    # Test multiple mutations to ensure variety
    samples = 500
    mutations = [mutate(setup.test_indiv_small) for _ in 1:samples]
    # Check that we get different mutations (very unlikely to get all identical)
    diff_mutes = length(unique(map(m -> m.ops, mutations)))
    @test diff_mutes > 1
    # should get some identical though (small sample)
    @test diff_mutes < (samples * .9)
    #print("Different mutes for ", samples," samples: ", diff_mutes,"\n")

end