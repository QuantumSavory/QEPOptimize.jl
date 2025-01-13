to run an example that simply evaluates the performance of a hardcoded circuit:

```
julia --project=. -tauto,auto example/evaluate.jl
```

to run an actual optimization

```
julia --project=. -tauto,auto example/optimize.jl
```