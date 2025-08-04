using TestItems

@testitem "cleanup_two_measurements"  begin
    using BPGates
    using QEPOptimize:cleanup_two_measurements!

    a = BellMeasure(1,5)
    b = CNOTPerm(1,2,3,4)
    c = BellMeasure(3,1)
    d = BellMeasure(3,1)

    o = [a,b,c,d]

    @test cleanup_two_measurements!(o,5) == [a,b,c]

    o = [a,b,c,d,d,d,d,d,d,d]
    @test cleanup_two_measurements!(o,5) == [a,b,c]

    o = [a,d]
    @test cleanup_two_measurements!(o,5) == [a,d]

    o = [b, b]
    @test cleanup_two_measurements!(o,5) == [b,b]

    # difficult case where the 'ops' are not back to back in the array, but they are in the actual circuit

    a = BellMeasure(1,2)
    b = CNOTPerm(1,2,3,4)
    c = BellMeasure(3,2)
    o = [a,b,c]
    @test cleanup_two_measurements!(o,5) == [a,b]

end

@testitem "cleanup_untargetted_pairs" begin
    using BPGates
    using QEPOptimize:get_used_qubits,cleanup_untargetted_pairs!

    a = BellMeasure(1,5)
    b = CNOTPerm(1,2,3,4)
    c = BellMeasure(3,1)
    d = BellMeasure(3,1)

    o = [a,b,c,d]

    @test get_used_qubits(o) == Set([1,3,4,5])

    # consider circuit on 3 pairs, with 1 purified 
    e = CNOTPerm(1,2,1,2)
    d = BellMeasure(1,2)

    # not using pair 3
    o = [e,d]
    @test get_used_qubits(o) == Set([1,2])

    cleanup_untargetted_pairs!(o,3,1)

    # Should have added a CNOT and measure on last pair
    @test get_used_qubits(o) == Set([1,2,3])
    # test that the correct gates were added (CNOTPerm with idx2 = 3, idx1 = 1, and BellMeasure after with sidx=3)

    @test o[1] isa CNOTPerm
    @test o[1].idx2 == 3 
    @test o[1].idx1 == 1
    @test o[2] isa BellMeasure 
    @test o[2].sidx == 3
        

    # Now consider if the purification qubit is not used 
    f = CNOTPerm(1,2,2,3)
    g = BellMeasure(1,3)
    o = [f,g]
    @test get_used_qubits(o) == Set([2,3])

    cleanup_untargetted_pairs!(o,3,1)

    # Should have added a CNOT with control qubit being the purified pair
    @test get_used_qubits(o) == Set([1,2,3])
    # test that the correct gate was added (cnot with idx1 = 1)

    @test o[1] isa CNOTPerm
    @test o[1].idx1 == 1
    

    # Make sure it works with empty circ 
    o = []
    cleanup_untargetted_pairs!(o,3,1)
    @test get_used_qubits(o) == Set([1,2,3])

    # does nothing on failure config 
    o = []
    cleanup_untargetted_pairs!(o,3,3)
    @test get_used_qubits(o) == Set()
end

@testitem "cleanup_measurements_on_top_qubits!"  begin
    using BPGates
    using QEPOptimize:cleanup_measurements_on_top_qubits!

    a = BellMeasure(1,5)
    b = CNOTPerm(1,2,3,4)
    c = BellMeasure(3,1)
    d = BellMeasure(3,1)

    o = [a,b,c,d]

    @test cleanup_measurements_on_top_qubits!(o,1) == [a,b]

    o = [a,b,c,d,d,d,d,d,d,d]
    @test cleanup_measurements_on_top_qubits!(o,1) == [a,b]

    o = [a,b,c,d]
    @test cleanup_measurements_on_top_qubits!(o,5) == [b]
end

@testitem "cleanup_nonmeasurement_last_steps!"  begin
    using BPGates
    using QEPOptimize:cleanup_nonmeasurement_last_steps!
    using QuantumClifford:AbstractMeasurement

    # 3 pairs, 1 pure, no last nonmeasurements
    a = CNOTPerm(1,2,3,2)
    b = BellMeasure(1,3)
    c = BellMeasure(3,2)

    o = [a,b,c]

    @test cleanup_nonmeasurement_last_steps!(o,3,1) == [a,b,c]

    # now let's say there are 4 pairs, but we did not even use the 4th -> no change should happen 

    @test cleanup_nonmeasurement_last_steps!(o,4,1) == [a,b,c]

    # now we add a cnot at the end from 3 to 4 which means 3 and 4 should have added measures 
    d = CNOTPerm(1,2,3,4) 
    o = [a,b,c,d]
    cleanup_nonmeasurement_last_steps!(o,4,1)
    @test length(o) == 6
    # does not matter the order, but the last two ops should both be measures on 3 and 4 
    @test typeof(o[5]) <: AbstractMeasurement
    @test typeof(o[6]) <: AbstractMeasurement
    @test (o[5].sidx == 3 && o[6].sidx == 4) || (o[5].sidx == 4 && o[6].sidx == 3)

    # now we add a cnot at the end from 1 to 4 which means 4 should have added measure 
    d = CNOTPerm(1,2,1,4) 
    o = [a,b,c,d]
    cleanup_nonmeasurement_last_steps!(o,4,1)
    @test length(o) == 5
    @test typeof(o[5]) <: AbstractMeasurement
    @test o[5].sidx == 4

    # and test that it still runs even though the check can't be run because of a bad config 
    o = [a,b,c,d]
    cleanup_nonmeasurement_last_steps!(o,4,4)
    @test length(o) == 4
end

@testitem "canonicalization is applied for optimization - simple" begin
    # run a simulation 
    using QEPOptimize
    using QEPOptimize: initialize_pop!, step!, cleanup_untargetted_pairs!,cleanup_measurements_on_top_qubits!,cleanup_nonmeasurement_last_steps!,cleanup_two_measurements!

    pop = Population()
    num_pairs = 4;num_purified = 1;steps = 5;pop_size = 20

    initialize_pop!(pop; number_registers=num_pairs,purified_pairs=num_purified,pop_size=pop_size)

    _, fitness_history, transition_counts_matrix, transition_counts_keys = multiple_steps_with_history!(pop, steps; number_registers=num_pairs,purified_pairs=num_purified,pop_size=pop_size)

    # all of the circuits should already be canonicalized
    for indiv in pop.individuals
        old_ops = copy(indiv.ops)
        @test cleanup_two_measurements!(indiv.ops,num_pairs) == old_ops
        @test cleanup_untargetted_pairs!(indiv.ops,num_pairs,num_purified) == old_ops
        @test cleanup_measurements_on_top_qubits!(indiv.ops,num_purified) == old_ops
        @test cleanup_nonmeasurement_last_steps!(indiv.ops,num_pairs,num_purified) == old_ops
    end
end