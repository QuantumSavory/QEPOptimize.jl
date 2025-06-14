### A Pluto.jl notebook ###
# v0.20.10

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
# ╠═╡ show_logs = false
begin
	using Pkg
	Pkg.pkg"add Revise, CairoMakie, PlutoUI, Quantikz, BPGates, QuantumClifford, ProgressLogging"
	#using Revise
	#Pkg.add(url="https://github.com/QuantumSavory/QEPOptimize.jl.git")
	Pkg.develop(path="../")
	using CairoMakie
	using PlutoUI
	using Quantikz
	using QEPOptimize
	using ProgressLogging
	using QEPOptimize: initialize_pop!, step!, NetworkFidelity, to_qasm, to_stabilizer
	using BPGates
	using BPGates: PauliNoise, BellMeasure, CNOTPerm
	using QuantumClifford: SparseGate, sCNOT, affectedqubits, BellMeasurement, Reset, sMX, sMZ, sMY
end

# ╔═╡ 8fc5cb18-70cc-4846-a62b-4cda69df12b0
md"
# QEPOptimize.jl
Entanglement Purification Circuit Generator
"

# ╔═╡ 6419143d-dc3a-47f0-8791-004e57b911c1
md"""
## Quantum Circuit Parameters

* Number of registers: $(@bind number_registers PlutoUI.Slider(2:6, default=4, show_value=true))

* Purified Pairs: $(@bind purified_pairs PlutoUI.Slider(1:5, default=1, show_value=true))

* Maximum Operations: $(@bind max_ops PlutoUI.Slider(10:5:30, default=15, show_value=true))

### Error Parameters

* Network fidelity: $(@bind network_fidelity PlutoUI.Slider(0.:0.002:0.3, default=0.1, show_value=true))

* Gate error X: $(@bind paulix PlutoUI.Slider(0.:0.002:0.1, default=0.01, show_value=true))

* Gate error Y: $(@bind pauliy PlutoUI.Slider(0.:0.002:0.1, default=0.01, show_value=true))

* Gate error Z: $(@bind pauliz PlutoUI.Slider(0.:0.002:0.1, default=0.01, show_value=true))

## Simulation Parameters

* Number of Simulations: $(@bind num_simulations PlutoUI.Slider(100:100:10000, default=1000, show_value=true))

* Population Size: $(@bind pop_size PlutoUI.Slider(10:10:100, default=20, show_value=true))

* Initial Operations: $(@bind start_ops PlutoUI.Slider(5:20, default=10, show_value=true))
* Initial Population Size: $(@bind start_pop_size PlutoUI.Slider(100:100:2000, default=1000, show_value=true))

### Evolution Parameters

* Number of Evolution Steps: $(@bind evolution_steps PlutoUI.Slider(50:150, default=80, show_value=true))

* New Mutants: $(@bind new_mutants PlutoUI.Slider(5:5:30, default=10, show_value=true))

* Drop Probability: $(@bind p_drop PlutoUI.Slider(0.0:0.05:0.5, default=0.1, show_value=true))

* Mutation Probability: $(@bind p_mutate PlutoUI.Slider(0.0:0.05:0.5, default=0.1, show_value=true))

* Gain Probability: $(@bind p_gain PlutoUI.Slider(0.0:0.05:0.5, default=0.1, show_value=true))

"""

# ╔═╡ 7419143d-dc3a-47f0-8791-004e57b911c2
begin
	config = (;
	    num_simulations=num_simulations,
	    number_registers=number_registers,
	    purified_pairs=purified_pairs,
	    code_distance=1, # Not needed to change for now
	    pop_size=pop_size,
	    noises=[NetworkFidelity(network_fidelity), PauliNoise(paulix, pauliy, pauliz)],
	)
	
	init_config = (;
	    start_ops=start_ops,
	    start_pop_size=start_pop_size,
	    config...
	)
	
	step_config = (;
	    max_ops=max_ops,
	    new_mutants=new_mutants,
	    p_drop=p_drop,
	    p_mutate=p_mutate,
	    p_gain=p_gain,
	    config...
	)
	
	pop = Population()
	
	initialize_pop!(pop; init_config...) 
end;


# ╔═╡ c09c7bb8-1d08-45da-81ca-0cf1d1985b91
_, fitness_history, transition_counts_matrix, transition_counts_keys = multiple_steps_with_history!(pop, evolution_steps; step_config...); 

# ╔═╡ 451be68d-b0bb-4b1b-b7fa-5c39618f95de
md"
## Simulation Results and Fidelity over generations
"

# ╔═╡ 988e9e99-cf93-46a3-be59-11c11e316b07
plot_fitness_history(
	fitness_history,
	transition_counts_matrix,
	transition_counts_keys
)

# ╔═╡ 1b6a9400-9d3b-42f1-a83f-c16f8134cb93
md"
## Best Circuit
"


# ╔═╡ e876ddcf-d2c9-401e-af83-368fbd5ba593
begin
	best_circuit = pop.individuals[1]
	Quantikz.QuantikzOp.(best_circuit.ops)
end

# ╔═╡ 4ab68db5-70cd-45e1-90bb-9fbb2830a3e4
md"
Fidelity results for this circuit

* Fidelity Simulations:  $(@bind fidelity_num_simulations PlutoUI.Slider(10000:10000:200000, default=100000, show_value=true))
"

# ╔═╡ 81aa21b4-50f0-4695-a9d0-fd998b0c0cc1
plot_circuit_analysis(
	best_circuit;
	num_simulations=fidelity_num_simulations,
	config.number_registers,
	config.purified_pairs,
    noise_sets=[[PauliNoise(paulix, pauliy, pauliz)],[]],
    noise_set_labels=["with gate noise", "no gate noise"]
)

# ╔═╡ 55dec933-ace8-4a90-bfc3-3dcd9e23a4cc
md"
Operations on this circuit:
"

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


# ╔═╡ d5fc459c-2c34-4f2f-a7ab-bce5e196ce1e
begin
	qasm = to_qasm(best_circuit.ops,number_registers,purified_pairs;comments=true)
	print(qasm)
end

# ╔═╡ 1a21ace6-b448-49d8-b544-97e8c0194723
s_out = to_stabilizer(best_circuit.ops,number_registers)

# ╔═╡ eac9bbe6-29c4-402b-84af-c88fb826a782
print(s_out)

# ╔═╡ 7674a28f-7bf2-426a-817f-68f83e6425d7


# ╔═╡ Cell order:
# ╟─8fc5cb18-70cc-4846-a62b-4cda69df12b0
# ╠═353e15de-0a9b-4107-a265-28953e1deee2
# ╟─6419143d-dc3a-47f0-8791-004e57b911c1
# ╟─7419143d-dc3a-47f0-8791-004e57b911c2
# ╠═c09c7bb8-1d08-45da-81ca-0cf1d1985b91
# ╟─451be68d-b0bb-4b1b-b7fa-5c39618f95de
# ╠═988e9e99-cf93-46a3-be59-11c11e316b07
# ╟─1b6a9400-9d3b-42f1-a83f-c16f8134cb93
# ╠═e876ddcf-d2c9-401e-af83-368fbd5ba593
# ╠═4ab68db5-70cd-45e1-90bb-9fbb2830a3e4
# ╠═81aa21b4-50f0-4695-a9d0-fd998b0c0cc1
# ╠═55dec933-ace8-4a90-bfc3-3dcd9e23a4cc
# ╠═e19cb382-99ae-4629-8242-83827c9e3631
# ╠═49894406-6dfe-4aeb-8193-e31731bfab65
# ╠═e99563b8-2827-4e0a-b7c5-e812fef7c6c5
# ╠═c434086a-d9e3-436b-91ad-a7ddef56622d
# ╠═d5fc459c-2c34-4f2f-a7ab-bce5e196ce1e
# ╠═1a21ace6-b448-49d8-b544-97e8c0194723
# ╠═eac9bbe6-29c4-402b-84af-c88fb826a782
# ╠═7674a28f-7bf2-426a-817f-68f83e6425d7
