# News

## v0.1.1 - 2025-07-17

- Pluto reactivity. Running the pluto notebook is a lot more enjoyable from the user side. Configuration changes are buffered by a PlutoUI.confirm query, meaning many variables can be changed before any simulation is run. Restarting/running the simulation also are now buffered with buttons.
- Pluto output select. The user can select which outputs to see

Things to note:
- To make the reactivity improvements on pluto, some hacks are used (see the second code block, and the config block). This increases load time by a small amount.


## v0.1.0 - 2025-06-13

- Simulation capabilities for commond gates.
- Simulation capabilities for Pauli noise.
- Simple genetic algorithm for optimization of circuits
- Simple Pluto UI example

