# QEPOptimize.jl - Quantum Entanglement Purification Optimizer 
The purpose of this package is to generate quantum circuits that create error-corrected entanglement, based on [this paper](https://arxiv.org/abs/2307.06354). The optimizer uses a genetic algorithm that edits circuits and evaluates their performance with [Monte Carlo trajectories](https://qc.quantumsavory.org/stable/API/#QuantumClifford.mctrajectory!-Tuple%7BAny,%20Any%7D). QEPOptimize.jl is built on top of previous work done by @Krastanov and @YipiaoWu, and relies on [QuantumSavory/BPGates.jl](https://github.com/QuantumSavory/BPGates.jl) and [QuantumSavory/QuantumClifford.jl](https://github.com/QuantumSavory/QuantumClifford.jl), also using visuals from [QuantumSavory/Quantikz.jl](https://github.com/QuantumSavory/Quantikz.jl).
## Pluto UI

To launch the Pluto notebook interface for configuring and running an optimization, run the following script in your terminal:

```sh
julia --project=. -tauto,auto -e 'using Pluto; Pluto.run(notebook="example/pluto.jl")'
```

If you haven't installed Pluto yet, add it first:

```julia
] add Pluto
```
## Run without pluto
To run an example that simply evaluates the performance of a hardcoded circuit:

```
julia --project=. -tauto,auto example/evaluate.jl
```

To run an actual optimization

```
julia --project=. -tauto,auto example/optimize.jl
```