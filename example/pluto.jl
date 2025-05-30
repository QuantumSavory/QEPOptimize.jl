### A Pluto.jl notebook ###
# v0.20.8

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 353e15de-0a9b-4107-a265-28953e1deee2
begin
	using Pkg
	Pkg.pkg"add Revise, CairoMakie, PlutoUI, Quantikz, BPGates, QuantumClifford"
	#using Revise
	#Pkg.add(url="https://github.com/QuantumSavory/QEPOptimize.jl.git")
	Pkg.develop(path="../")
	using CairoMakie
	using PlutoUI
	using Quantikz
	using QEPOptimize
	using QEPOptimize: initialize_pop!, step!, NetworkFidelity
	using BPGates
	using BPGates: PauliNoise, BellMeasure, CNOTPerm
	using QuantumClifford: SparseGate, sCNOT, affectedqubits, BellMeasurement, Reset, sMX, sMZ, sMY
end

# ╔═╡ be98a0ab-b122-4a88-a02d-d57dad7cbc13
md"""
number of registers: $(@bind number_registers PlutoUI.Slider(2:6, default=4, show_value=true))

gate error X: $(@bind paulix PlutoUI.Slider(0.:0.002:0.1, default=0.01, show_value=true))

gate error Y: $(@bind pauliy PlutoUI.Slider(0.:0.002:0.1, default=0.01, show_value=true))

gate error Z: $(@bind pauliz PlutoUI.Slider(0.:0.002:0.1, default=0.01, show_value=true))

network fidelity: $(@bind network_fidelity PlutoUI.Slider(0.:0.002:0.3, default=0.1, show_value=true))
"""

# ╔═╡ 6419143d-dc3a-47f0-8791-004e57b911c1
begin
	config = (;
	    num_simulations=1000,
	    number_registers=number_registers,
	    purified_pairs=1,
	    code_distance=1,
	    pop_size = 20,
	    noises=[NetworkFidelity(network_fidelity), PauliNoise(paulix, pauliy, pauliz)],
	)
	
	init_config = (;
	    start_ops = 10,
	    start_pop_size = 1000,
	    config...
	)
	
	step_config = (;
	    max_ops = 15, # do 5 for something faster
	    new_mutants = 10,
	    p_drop = 0.1,
	    p_mutate = 0.1,
	    p_gain = 0.1,
	    config...
	)
	
	pop = Population()
	
	initialize_pop!(pop; init_config...)
end;

# ╔═╡ c09c7bb8-1d08-45da-81ca-0cf1d1985b91
_, fitness_history, transition_counts_matrix, transition_counts_keys = multiple_steps_with_history!(pop, 80; step_config...);

# ╔═╡ 988e9e99-cf93-46a3-be59-11c11e316b07
plot_fitness_history(
	fitness_history,
	transition_counts_matrix,
	transition_counts_keys
)

# ╔═╡ e876ddcf-d2c9-401e-af83-368fbd5ba593
begin
	best_circuit = pop.individuals[1]
	Quantikz.QuantikzOp.(best_circuit.ops)
end

# ╔═╡ 81aa21b4-50f0-4695-a9d0-fd998b0c0cc1
plot_circuit_analysis(
	best_circuit;
	num_simulations=100000,
	config.number_registers,
	config.purified_pairs,
    noise_sets=[[PauliNoise(paulix, pauliy, pauliz)],[]],
    noise_set_labels=["with gate noise", "no gate noise"]
)

# ╔═╡ e19cb382-99ae-4629-8242-83827c9e3631
begin
	for g in best_circuit.ops
		println(g)
	end
end

# ╔═╡ 49894406-6dfe-4aeb-8193-e31731bfab65
@doc BPGates.CNOTPerm

# ╔═╡ e99563b8-2827-4e0a-b7c5-e812fef7c6c5
@doc BPGates.BellMeasure

# ╔═╡ c434086a-d9e3-436b-91ad-a7ddef56622d
begin
	for g in best_circuit.ops
		for qcg in BPGates.toQCcircuit(g)
			if isa(qcg, SparseGate)
				println(qcg.cliff)
				println("on qubit $(qcg.indices...)")
			elseif isa(qcg, sCNOT)
				println("CNOT on qubits $(affectedqubits(qcg))")
			elseif isa(qcg, BellMeasurement)
				print("measure ")
				for m in qcg.measurements
				    print(Dict(sMX=>:X, sMY=>:Y, sMZ=>:Z)[typeof(m)])
					print(m.qubit)
					print(" ")
				end
				println("of parity $(qcg.parity)")
			elseif isa(qcg, Reset)
				println("new raw Bell pair")
			else
				println(qcg)
			end
			println()
		end
	end
end

# ╔═╡ Cell order:
# ╠═353e15de-0a9b-4107-a265-28953e1deee2
# ╠═be98a0ab-b122-4a88-a02d-d57dad7cbc13
# ╠═6419143d-dc3a-47f0-8791-004e57b911c1
# ╟─c09c7bb8-1d08-45da-81ca-0cf1d1985b91
# ╠═988e9e99-cf93-46a3-be59-11c11e316b07
# ╟─e876ddcf-d2c9-401e-af83-368fbd5ba593
# ╠═81aa21b4-50f0-4695-a9d0-fd998b0c0cc1
# ╠═e19cb382-99ae-4629-8242-83827c9e3631
# ╠═49894406-6dfe-4aeb-8193-e31731bfab65
# ╠═e99563b8-2827-4e0a-b7c5-e812fef7c6c5
# ╠═c434086a-d9e3-436b-91ad-a7ddef56622d
