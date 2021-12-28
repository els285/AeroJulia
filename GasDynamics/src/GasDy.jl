module GasDynamics

include("BaseFluid.jl")
using .BaseFluid

include("Isentropic.jl")
using .Isentropic

include("Shocks.jl")
using .Shocks

include("./PMExpansion.jl")
using .PMExpansion

include("./ConicalFlow.jl")
using .ConicalFlow


end