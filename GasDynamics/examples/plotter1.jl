# Example PLot for Julia
include("GasDy.jl")
using .GasDynamics

# include("PlotRatios.jl")
# using .GasDynamic_Plots


f1 = GasDynamics.Fluid("air",1.4,287)


# min_beta = GasDynamics.Shocks.ObliqueShock.Minimum_ShockAngle(1.6,1.4)
# println(min_beta*180/pi)
GasDynamics.Shocks.Plot_ThetaBetaM(f1)
# GasDynamics.Shocks.Plot_ObliqueRatios_beta(f1,3)
# println(GasDynamics.Shocks.NormalRatios)


# GasDynamic_Plots.Make_Ratio_Plot(f1,"Isentropic",M2=0.5,M_range=(1,5))
# include("Shocks.jl")
# using .Shocks

# println(f1)