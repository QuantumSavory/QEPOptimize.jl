### A Pluto.jl notebook ###
# v0.20.13

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
	using PlutoUI:confirm, Slider
	using Quantikz
	using QEPOptimize
	using QEPOptimize: initialize_pop!, step!, NetworkFidelity, Population, EVOLUTION_METRICS
	using BPGates
	using BPGates: PauliNoise, BellMeasure, CNOTPerm
	using QuantumClifford: SparseGate, sCNOT, affectedqubits, BellMeasurement, Reset, sMX, sMZ, sMY
end

# ╔═╡ 8fc5cb18-70cc-4846-a62b-4cda69df12b0
md"
# QEPOptimize.jl
Entanglement Purification Circuit Generator
"

# ╔═╡ 610c0135-3a8e-4676-aa5f-9ca76546dd98
begin
	# to help dealing with variables that should not be changed ever
	
	# Define the population (on its own so it is not regenerated)
	const pop = Ref(Population())

	# These are here to allow re-definition of the configs' values, without triggering cells (only change the derefrenced value) 
	# the types are nasty... I know..
	
	const config = Ref{@NamedTuple{num_simulations::Int64, number_registers::Int64, purified_pairs::Int64, code_distance::Int64, pop_size::Int64, noises::Vector{Any}}}()
	
	const init_config = Ref{@NamedTuple{start_ops::Int64, start_pop_size::Int64, num_simulations::Int64, number_registers::Int64, purified_pairs::Int64, code_distance::Int64, pop_size::Int64, noises::Vector{Any}}}()

	const step_config = Ref{@NamedTuple{max_ops::Int64, new_mutants::Int64, p_drop::Float64, p_mutate::Float64, p_gain::Float64, evolution_metric::Symbol, max_performance_calcs::Int64, max_simulations::Int64, num_simulations::Int64, number_registers::Int64, purified_pairs::Int64, code_distance::Int64, pop_size::Int64, noises::Vector{Any}}}()

	const evolution_steps_ref = Ref{Int64}()

	nothing;
end

# ╔═╡ 6419143d-dc3a-47f0-8791-004e57b911c1
@bind c confirm(
	PlutoUI.combine() do Child
md"""
## Quantum Circuit Parameters

* Number of registers: $(Child("number_registers", PlutoUI.Slider(2:6, default=4, show_value=true)))

* Purified Pairs: $(Child("purified_pairs", PlutoUI.Slider(1:5, default=1, show_value=true)))

* Maximum Operations: $(Child("max_ops", PlutoUI.Slider(10:5:30, default=15, show_value=true)))

* Code distance: $(Child("code_distance", PlutoUI.Slider(1:1:6, default=1, show_value=true)))

### Error Parameters

* Network fidelity: $(Child( "network_fidelity", PlutoUI.Slider(0.:0.002:0.3, default=0.1, show_value=true)))

* Gate error X: $(Child("paulix",PlutoUI.Slider(0.:0.002:0.1, default=0.01, show_value=true)))

* Gate error Y: $(Child("pauliy", PlutoUI.Slider(0.:0.002:0.1, default=0.01, show_value=true)))

* Gate error Z: $(Child( "pauliz",PlutoUI.Slider(0.:0.002:0.1, default=0.01, show_value=true)))

## Simulation Parameters

* Number of Simulations: $(Child("num_simulations", PlutoUI.Slider(100:100:5000, default=1000, show_value=true)))

* Max number of Simulations: $(Child("max_simulations", PlutoUI.Slider(100:100:10000, default=5000, show_value=true)))

* Max performance calculations per circuit: $(Child("max_perf_calcs", PlutoUI.Slider(1:1:50, default=10, show_value=true)))

* Population Size: $( Child("pop_size", PlutoUI.Slider(10:10:100, default=20, show_value=true)))

* Initial Operations: $(Child( "start_ops", PlutoUI.Slider(5:20, default=10, show_value=true)))
* Initial Population Size: $(Child( "start_pop_size", PlutoUI.Slider(100:100:2000, default=1000, show_value=true)))

### Evolution Parameters

* Evolution metric: $(Child("evolution_metric", PlutoUI.Select([:logical_qubit_fidelity => "Logical qubit fidelity",:purified_pairs_fidelity => "Purified pairs fidelity" ,:average_marginal_fidelity => "Average marginal fidelity"])))

* Number of Evolution Steps: $(Child("evolution_steps", PlutoUI.Slider(10:150, default=50, show_value=true)))

* New Mutants: $(Child( "new_mutants", PlutoUI.Slider(5:5:30, default=10, show_value=true)))

* Drop Probability: $(Child("p_drop", PlutoUI.Slider(0.0:0.05:0.5, default=0.1, show_value=true)))

* Mutation Probability: $(Child("p_mutate", PlutoUI.Slider(0.0:0.05:0.5, default=0.1, show_value=true)))

* Gain Probability: $(Child("p_gain", PlutoUI.Slider(0.0:0.05:0.5, default=0.1, show_value=true)))

"""
	end; label="Update Config"
)

# ╔═╡ a892e297-7223-4d1a-b772-5f4ca5c64339
@bind restart_population_trigger PlutoUI.CounterButton("Restart Population")

# ╔═╡ dad1728c-c341-44cc-88e6-d26ca1815a30
@bind run_simulation_trigger PlutoUI.CounterButton("Run simulation")

# ╔═╡ cef70317-fc58-42b3-987b-a454064f0113
begin
	config[] = (
			num_simulations=c.num_simulations,
		    number_registers=c.number_registers,
		    purified_pairs=c.purified_pairs,
		    code_distance=c.code_distance, # For logical qubit fidelity
		    pop_size=c.pop_size,
			noises=[NetworkFidelity(c.network_fidelity), PauliNoise(c.paulix, c.pauliy, c.pauliz)],
		)
	
	init_config[] = (
		start_ops=c.start_ops,
		start_pop_size=c.start_pop_size,
		config[]...
	)

	step_config[] = (;
		max_ops=c.max_ops,
		new_mutants=c.new_mutants,
		p_drop=c.p_drop,
		p_mutate=c.p_mutate,
		p_gain=c.p_gain,
		evolution_metric=c.evolution_metric,
		max_performance_calcs=c.max_perf_calcs,
		max_simulations=c.max_simulations,
		config[]...
	)

	evolution_steps_ref[] = c.evolution_steps
end;

# ╔═╡ 3d17bc74-fa91-410c-b060-b15eae7a564b
begin
	# Re-generate population
	restart_population_trigger
	
	initialize_pop!(pop[]; init_config[]...); 
	md"Regenerating Population..."
end

# ╔═╡ c09c7bb8-1d08-45da-81ca-0cf1d1985b91
begin 
	# Run simulation
	run_simulation_trigger
	try
		# check if values have been set
		pop[], evolution_steps_ref[], step_config[];
	catch
		throw("Config not set")
	end
	_, fitness_history, transition_counts_matrix, transition_counts_keys = multiple_steps_with_history!(pop[], evolution_steps_ref[]; step_config[]...); 
	
	md"""Optimizing..."""
end

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
	best_circuit = pop[].individuals[1]
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
	config[].number_registers,
	config[].purified_pairs,
    noise_sets=[[PauliNoise(c.paulix, c.pauliy, c.pauliz)],[]],
    noise_set_labels=["with gate noise", "no gate noise"]
)


# ╔═╡ 402e6b8b-7c13-4ab0-9d40-1b764dba1691
md"
Choose how to display operations on this circuit: 
"

# ╔═╡ 9a1b169a-36f5-4d2e-ad63-a7e184abde66
begin
	@bind enabled_output PlutoUI.MultiCheckBox([:print => "Print",:stab_desc => "Stabilizer description"];default=[:print,:stab_desc])
end

# ╔═╡ 23123ce9-58b0-4eb7-8d39-fd56499b3ed2
:print in enabled_output ? md"#### Print" : nothing

# ╔═╡ e19cb382-99ae-4629-8242-83827c9e3631
begin
	if :print in enabled_output
		for g in best_circuit.ops println(g) end
		md"""
		BPGates.CNOTPerm:
		$(@doc BPGates.CNOTPerm)
		
		BPGates.BellMeasure:
		$(@doc BPGates.BellMeasure)

		Operations on the best circuit:
		"""
	else
	end
end

# ╔═╡ b2f3c15c-1b03-4c4e-9c68-539f96ebd4cb
:stab_desc in enabled_output ? md"#### Stabilizer Description" : nothing

# ╔═╡ c434086a-d9e3-436b-91ad-a7ddef56622d
begin
	if :stab_desc in enabled_output
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
end

# ╔═╡ Cell order:
# ╟─8fc5cb18-70cc-4846-a62b-4cda69df12b0
# ╠═353e15de-0a9b-4107-a265-28953e1deee2
# ╠═610c0135-3a8e-4676-aa5f-9ca76546dd98
# ╟─6419143d-dc3a-47f0-8791-004e57b911c1
# ╟─a892e297-7223-4d1a-b772-5f4ca5c64339
# ╟─dad1728c-c341-44cc-88e6-d26ca1815a30
# ╠═cef70317-fc58-42b3-987b-a454064f0113
# ╟─3d17bc74-fa91-410c-b060-b15eae7a564b
# ╠═c09c7bb8-1d08-45da-81ca-0cf1d1985b91
# ╟─451be68d-b0bb-4b1b-b7fa-5c39618f95de
# ╟─988e9e99-cf93-46a3-be59-11c11e316b07
# ╟─1b6a9400-9d3b-42f1-a83f-c16f8134cb93
# ╟─e876ddcf-d2c9-401e-af83-368fbd5ba593
# ╟─4ab68db5-70cd-45e1-90bb-9fbb2830a3e4
# ╟─81aa21b4-50f0-4695-a9d0-fd998b0c0cc1
# ╟─402e6b8b-7c13-4ab0-9d40-1b764dba1691
# ╟─9a1b169a-36f5-4d2e-ad63-a7e184abde66
# ╟─23123ce9-58b0-4eb7-8d39-fd56499b3ed2
# ╟─e19cb382-99ae-4629-8242-83827c9e3631
# ╟─b2f3c15c-1b03-4c4e-9c68-539f96ebd4cb
# ╟─c434086a-d9e3-436b-91ad-a7ddef56622d
