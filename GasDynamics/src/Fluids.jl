# Fluids


include("BaseFluid.jl")
using .BaseFluid


Air = BaseFluid.MakeFluid("Standard dry air, 293K",1.4,287)