include("./../src/GasDy.jl")
using .GasDynamics

f1 = GasDynamics.Fluid("air",1.4,287)

flow = GasDynamics.MakeFlow(f1,17;Density=1.125,Pressure=1e5)

# downstream_normal = GasDynamics.ApplyNormalShock(flow)

# println(downstream_normal.Pressure/(downstream_normal.Temperature*downstream_normal.Density))

# println(GasDynamics.Isentropic.Generate_IsentropicRatios(3,2,1.4))
# println(GasDynamics.ApplyExpansionFan(flow,0.5).Pressure)

downstream_oblique = GasDynamics.ApplyObliqueShock(flow;beta=0.3)

println(downstream_oblique)

# println(downstream_oblique.Pressure/(downstream_oblique.Temperature*downstream_oblique.Density))
# 
# x = GasDynamics.ObliqueRatios(2.5,1.4,1)



# println(flow)