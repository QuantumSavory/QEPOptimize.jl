# TODO (low priority) this would be a great place to use an Enum or, even better, an algebraic data type (ADT)
const HISTORIES = [:manual, :survivor, :random, :child, :drop, :gain, :swap, :opsmutate]

"A convenient structure to store various purification performance metrics."
struct Performance
    "a list of probabilities as such `[probability for no errors, probability for one Bell pair to be erroneous, probability for two Bell pairs to be erroneous, ..., probability for all Bell pairs to be erroneous]`"
    error_probabilities::Vector{Float64} # TODO (low priority) this being a parameterized Tuple would make it statically inlineable
    "probability for no errors"
    purified_pairs_fidelity::Float64
    "the probability for no logical error if the purified pairs were used for teleportation of a code whose distance is a configuration option given to the performance evaluator/optimizer"
    logical_qubit_fidelity::Float64
    """the average marginal fidelity of all purified pairs (correlation of errors is ignored here);
    Fidelity of i-th Bell pair is `Fᵢ = ⟨A|Trᵢ(ρ)|A⟩` where `Trᵢ` is "partial trace for all subspaces except the i-th one".
    The average marginal fidelity is `mean([F₁, F₂, ... Fₖ])`."""
    average_marginal_fidelity::Float64
    "the proportion of runs of a given protocol that do not have detected errors (i.e. Alice and Bob do not measure any errors)"
    success_probability::Float64
end

Performance() = Performance(Float64[], 0.0, 0.0, 0.0, 0.0)

Base.copy(p::Performance) = Performance(copy(p.error_probabilities), p.purified_pairs_fidelity, p.logical_qubit_fidelity, p.average_marginal_fidelity, p.success_probability)


"An individual (circuit) in the population we are evolving"
mutable struct Individual
    "How did this individual come to be (conventionally a symbol from the `HISTORIES` list)"
    history::Symbol
    "A vector containing a sequence of quantum operations that make up the individual's circuit"
    ops::Vector{Any}
    "A variety of performance estimates"
    performance::Performance
    "Overall fitness, derived from `performance` and some configuration options given to the performance evaluator/optimizer"
    fitness::Float64
    # TODO (low priority) inner constructor that checks that `history` is a member of HISTORIES (will be unnecessary if we used an enum of an ADT)
end

Individual() = Individual(:manual)
Individual(history::Symbol) = Individual(history, [])
Individual(ops::Vector) = Individual(:manual, ops, Performance(), 0.0)
Individual(history::Symbol, ops) = Individual(history, ops, Performance(), 0.0)

Base.copy(i::Individual) = Individual(i.history, copy(i.ops), copy(i.performance), i.fitness)


mutable struct Population
    "All individuals in the population"
    individuals::Vector{Individual}
    "Keeps track of the selection history for different types of individuals (e.g., survivors, mutants)"
    selection_history::Dict{Symbol,Vector{Int64}}
end

Population() = Population([], Dict{Symbol, Vector{Int64}}())
