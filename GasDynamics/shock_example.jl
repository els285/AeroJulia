include("GasDy.jl")
using .GasDynamics

f1 = GasDynamics.Fluid("air",1.4,287)

flow = GasDynamics.MakeFlow(f1,2.5;Density=1.125,Pressure=1e5)

# downstream_normal = GasDynamics.ApplyNormalShock(flow)

# println(downstream_normal.Pressure/(downstream_normal.Temperature*downstream_normal.Density))

downstream_oblique = GasDynamics.ApplyObliqueShock(flow,0.3)

println(downstream_oblique.Pressure/(downstream_oblique.Temperature*downstream_oblique.Density))
# 
# x = GasDynamics.ObliqueRatios(2.5,1.4,1)



# println(flow)