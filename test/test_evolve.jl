using TestItems

@testitem "simulation can run if pop_size changes"  begin
    using QEPOptimize:multiple_steps_with_history!,Population,initialize_pop!

    pop_size = 10
    pop = Population()
    initialize_pop!(pop; pop_size=pop_size)
    @test size(pop.individuals)[1] == pop_size

    multiple_steps_with_history!(pop,5;pop_size=pop_size)
    @test size(pop.individuals)[1] == pop_size

    pop_size = 5
    multiple_steps_with_history!(pop,5;pop_size=pop_size)
    @test size(pop.individuals)[1] == pop_size

    pop_size = 20
    multiple_steps_with_history!(pop,5;pop_size=pop_size)
    @test size(pop.individuals)[1] == pop_size
end
