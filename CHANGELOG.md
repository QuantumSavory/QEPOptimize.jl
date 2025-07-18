# News

## v0.1.1 - 2025-07-17

- Qasm output. to_qasm function (added to pluto, with download button) converts any circuit to qasm counterpart. Testing needed, ideally with local qasm simulator.
- Stabilizer output. This ports the functionality from QuantumClifford, so this is a pretty small function. Also included in the pluto output.
- Pluto reactivity. Running the pluto notebook is a lot more enjoyable from the user side. Configuration changes are buffered by a PlutoUI.confirm query, meaning many variables can be changed before any simulation is run. Restarting/running the simulation also are now buffered with buttons.
- Pluto output select. The user can select which outputs to see, since there are now four total ways - the output can be quite a lot

Things to note:
- Testing is needed for the new functions and outputs
- Some improvements for qasm output and code style are currently waiting on https://github.com/QuantumSavory/BPGates.jl/pull/31 , 
- To make the reactivity improvements on pluto, some hacks are used (see the second code block, and the config block). This increases load time by a small amount.


## v0.1.0 - 2025-06-13

- Simulation capabilities for commond gates.
- Simulation capabilities for Pauli noise.
- Simple genetic algorithm for optimization of circuits
- Simple Pluto UI example

