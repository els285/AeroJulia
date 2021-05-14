# Load the GasDynamics module
include("GasDy.jl")
using .GasDynamics

# Generate a fluid
f1 = GasDynamics.Fluid("air",1.4,287)

# Define flow characteristics
flow = GasDynamics.MakeFlow(f1,2.5;Density=1.125,Pressure=1e5)

#Apply to normal shock wave
downstream_normal_shock = GasDynamics.ApplyNormalShock(flow)

# Apply to oblique shock wave
downstream_oblique_shock = GasDynamics.ApplyObliqueShock(flow,0.3)

# Apply to Prandtl-Meyer expansion fan
downstream_prandtlmeyer = GasDynamics.ApplyExpansionFan(flow,0.5)

