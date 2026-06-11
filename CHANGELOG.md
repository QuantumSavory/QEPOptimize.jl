# News

## v0.1.3 - 2025-08-22

- Generated circuits are 'cleaned up': canonicalized and edited if something is wrong. The first part of the optimization will run the 'unsafe' + 'safe' cleanups, and the second part only does the 'safe' ones. 
    - 'safe' cleanups are when there is no risk of making the circuit worse, such as removing excess measurements on the same pairs, or removing measurements on the purified pairs.
    - 'unsafe' cleanups risk making a circuit worse, but can also solve issues that the circuit may have, such as 'using' an unused pair by adding random a CNOTPerm and measurement on it, or adding a measurement at the end of a pair's ops if there are none. 
- The switch from unsafe+safe to only safe happens half way through (Likely to change)
- Adds Quantikz, BPGates 1.2.2 as deps

## v0.1.2 - 2025-08-12

- Pluto reactivity. Running the pluto notebook is a lot more enjoyable from the user side. Configuration changes are buffered by a PlutoUI.confirm query, meaning many variables can be changed before any simulation is run. Restarting/running the simulation also are now buffered with buttons.
- Pluto output select. The user can select which outputs to see

Things to note:
- To make the reactivity improvements on pluto, some hacks are used (see the second code block, and the config block). This increases load time by a small amount.

## v0.1.1 - 2025-07-29
- Performance evaluations carry over from previous runs. 
- Warning when simulation is throttled (fidelity = 1)

## v0.1.0 - 2025-06-13

- Simulation capabilities for commond gates.
- Simulation capabilities for Pauli noise.
- Simple genetic algorithm for optimization of circuits
- Simple Pluto UI example

