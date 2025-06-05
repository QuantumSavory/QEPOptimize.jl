function analyze_f_out_vs_f_in(
    circuit::Individual;
    num_simulations::Int=100000,
    number_registers::Int=2,
    purified_pairs::Int=1,
    noises=[PauliNoise(0.01/3, 0.01/3, 0.01/3)],
    f_ins = [0.01; 0.05:0.05:0.95; 0.99; 0.999]
)
    f_outs = Float64[]
    probs = Float64[]
    for f in f_ins
        p = calculate_performance!(circuit; num_simulations, noises=[NetworkFidelity(f),noises...], number_registers, purified_pairs)
        push!(f_outs, p.purified_pairs_fidelity)
        push!(probs, p.success_probability)
    end
    return f_ins, f_outs, probs
end

function plot_circuit_analysis(
    circuit::Individual;
    num_simulations::Int=100000,
    number_registers::Int=2,
    purified_pairs::Int=1,
    noise_sets=[[PauliNoise(0.01/3, 0.01/3, 0.01/3)],[]],
    noise_set_labels=[join(string.(noises), " ") for noises in noise_sets],
    f_ins = [0.01; 0.05:0.05:0.95; 0.99; 0.999]
)
    fig = Figure()
    axF = Axis(fig[1,1])
    axF.title = "F in vs F out"
    lines!(axF, [0,1], [0,1], color=:gray50)
    axP = Axis(fig[1,2])
    axP.title = "F in vs P"
    for (noises, label) in zip(noise_sets, noise_set_labels)
        f_ins, f_outs, probs = analyze_f_out_vs_f_in(circuit; num_simulations, number_registers, purified_pairs, noises, f_ins)
        lines!(axF, f_ins, f_outs)
        lines!(axP, f_ins, probs; label)
    end
    axislegend(axP, "Noise Model", position=:lt)
    return fig
end

function plot_fitness_history(
    fitness_history, transition_counts_matrix, transition_counts_keys
)
    fig = Figure(size=(600, 600))

    # heatmap of fitness history
    f_hist = fig[1,1]
    a_hist = Axis(f_hist[1,1], xlabel="Generation", ylabel="Individual")
    p_hist = plot!(a_hist, fitness_history[:,end:-1:begin], colorrange=(0.9, 1.0))
    c_hist = Colorbar(f_hist[1,2], p_hist, label="Fitness")

    # best and worst fitness over generations
    f_bestworst = fig[2,1]
    a_bestworst = Axis(f_bestworst[1,1], xlabel="Generation", ylabel="Fitness")
    lines!(a_bestworst, maximum(fitness_history, dims=2)[:], label="Best")
    lines!(a_bestworst, minimum(fitness_history, dims=2)[:], label="Worst")
    axislegend(a_bestworst, position=:rb)

    # type of circuit
    f_type = fig[3,1]
    a_type = Axis(f_type[1,1], xlabel="Generation", ylabel="Type")
    for i in 1:size(transition_counts_matrix, 2)
        lines!(a_type, transition_counts_matrix[:,i], label=string(transition_counts_keys[i]))
    end
    axislegend(a_type, position=:rb,labelsize=10,patchsize=(10,10))

    return fig
end
