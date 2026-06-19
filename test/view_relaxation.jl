using QEPOptimize: Individual,calculate_performance!
using BPGates: T1NoiseOp,T2NoiseOp

times = 0:1:2000
t1s = [100]
t1 = 100
λ₁_func(gate_time,t1) = 1-exp(-gate_time/(t1))
sims = 3000

function single_gate_data(tn,times,op_func=T1NoiseOp,lambda_func=λ₁_func;code_distance=1)
    logical_fidelities = Float64[]
    for time in times
        λ = lambda_func(time,tn)
        Tn_gate = [op_func(1, λ)]
        indiv_Tn = Individual(Tn_gate)
        perf_Tn = calculate_performance!(indiv_Tn; num_simulations=sims, noises=[], code_distance=code_distance)
        push!(logical_fidelities, perf_Tn.logical_qubit_fidelity)
    end
    return logical_fidelities
end

function many_gate_data(tn,dt,range, op_func=T1NoiseOp, lambda_func=λ₁_func;code_distance=1)
    logical_fidelities = Float64[]
    step = 0
    λ = lambda_func(dt,tn)
    total_steps = range / dt
    Tn_gate = []
    while step <= total_steps
        if step == 0
            push!(logical_fidelities, 1.0)
            step += 1
            continue
        else 
            append!(Tn_gate, [op_func(1, λ)]) 
            indiv_Tn = Individual(Tn_gate)
            perf_Tn = calculate_performance!(indiv_Tn; num_simulations=sims, noises=[], code_distance=code_distance)
            push!(logical_fidelities, perf_Tn.logical_qubit_fidelity)
            step += 1
        end
    end
    return logical_fidelities
end

function analytical_T1_logical_fidelity(t1,time)
    return 0.5 + 0.5*exp(-time/(3*t1))
end

using GLMakie
fig = Figure(resolution = (1200, 900))
ax = Axis(fig[1, 1]; xlabel = "Time (μs)", ylabel = "Logical Fidelity")

# Plot a line for evaluating one noise gate at a time
logical_fidelities = single_gate_data(t1, times)
lines!(ax, times, logical_fidelities; label = "single gate T1 = $(t1) μs")

# Plot lines for evaluating a circuit of many noise gates at a time
dt = 10
logical_fidelities_many = many_gate_data(t1, dt, maximum(times))
lines!(ax, 0:dt:maximum(times), logical_fidelities_many; label = "many gate T1 = $(t1) μs, dt=$(dt) μs")


dt = 50
logical_fidelities_many = many_gate_data(t1, dt, maximum(times), div(maximum(times), dt))
lines!(ax, 0:dt:maximum(times), logical_fidelities_many; label = "many gate T1 = $(t1) μs, dt=$(dt) μs")

analytical_fidelities = [analytical_T1_logical_fidelity(t1, time) for time in times]
lines!(ax, times, analytical_fidelities; linestyle = :dash, label = "analytical T1 = $(t1) μs")

ax.title = "Single-gate and multi-gate T1 Logical Fidelity"
axislegend(ax; position = :rt)
display(fig)




# T2 results

t1s = [100]
t2s = [150]

times = 0:1:2000
λ₂_func(gate_time,t1,t2) = 1 - exp(-gate_time/t2 + gate_time/(2t1))
sims = 3000


fig2 = Figure(resolution = (1200, 900))
ax2 = Axis(fig2[1, 1]; xlabel = "Time (μs)", ylabel = "Logical Fidelity")
for (t1,t2) in zip(t1s,t2s)
    lambda_func = (gate_time,tn) -> λ₂_func(gate_time,t1,t2)

    logical_fidelities = single_gate_data(t1, times, T2NoiseOp, lambda_func)
    lines!(ax2, times, logical_fidelities; label = "single gate T2 = $(t2) μs")

    dt = 10
    logical_fidelities_many = many_gate_data(t1, dt, maximum(times), div(maximum(times), dt),  T2NoiseOp   , lambda_func)
    lines!(ax2, 0:dt:maximum(times), logical_fidelities_many; label = "many gate T2 = $(t2) μs, dt=$(dt) μs")
   
    dt = 50
    logical_fidelities_many = many_gate_data(t1, dt, maximum(times), div(maximum(times), dt),  T2NoiseOp   , lambda_func)
    lines!(ax2, 0:dt:maximum(times), logical_fidelities_many; label = "many gate T2 = $(t2) μs, dt=$(dt) μs")

    # analytical_fidelities = [analytical_T1_logical_fidelity(t1, time) for time in times]
    # lines!(ax, times, analytical_fidelities; linestyle = :dash, label = "analytical T1 = $(t1) μs")


end
ax2.title = "Single-gate and multi-gate T2 Logical Fidelity, [T1=$(t1s[1]) μs]"
axislegend(ax2; position = :rt)
display(fig2)
