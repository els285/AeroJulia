# AeroJulia

A repository for aerospace-related scripts designed to help Julia users in aerospace-related work. Primarily focuses on gas dyanmics at the moment.


## BaseFluid

We use specific fluid structures to keep track of flow and fluid properties.

The fundamental object is a `Fluid`, which defines (as you'd guess) a fluid with a set of intrinsic mechanical and thermodynamic properties.

A second important object is a `Flow` which is the pairing of a `Fluid` with a particular kinematic and thermodynamic state (the MachNumber, velocity, pressure, temperature and density are all defined).
