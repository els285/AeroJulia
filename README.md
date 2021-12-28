# AeroJulia

A repository for aerospace-related scripts designed to help Julia users in aerospace-related work. 

Currently focuses solely on techniques for compressible flow problems.

# GasDynamics

## BaseFluid

We use specific fluid structures to keep track of flow and fluid properties.

The fundamental object is a `Fluid`, which defines (as you'd guess) a fluid with a set of intrinsic mechanical and thermodynamic properties. 

A second important object is a `Flow` which is the pairing of a `Fluid` with a particular kinematic and thermodynamic state (the MachNumber, velocity, pressure, temperature and density are all defined).

## Compressible Flows
The `GasDynamics` module defines a set of structures and functions useful for solving gas-dynamical and compressible flow problems. 

### Shocks
Capability to derive flow conditions behind **normal shock waves** and **oblique shock waves**, and apply these techniques to simple geometries in both internal and external flow contexts. 

### Expansions
Capability to derive flow conditions under **isentropic expansion** and under **Prandtl-Meyer expansion fans**, and apply these technqiues to simple geometries in both internal and external flow contexts.
