
module PMExpansion

using Roots

include("./Isentropic.jl")
using .Isentropic

include("./BaseFluid.jl")
using .BaseFluid

export ApplyExpansionFan , ApplyExpansionFan_Wrapped

"""
The Prandtl-Meyer expansion fan describes the flow of a quasi-one-directional compressible fluid as it expands isentropically around a corner

Here we define the Prandtl-Meyer function in radians

The package Roots is used to solve the Prandtl-Meyer function for Mach Number

"""


mutable struct PM_Wrapper
    """
    Wrapper for PM Expansion
    """
    Upstream 
    Downstream
    incoming_angle
    outgoing_angle 
end

function PrandtlMeyer(M,gamma)

    a = sqrt((gamma + 1) / (gamma - 1)) * atan(sqrt((gamma - 1) / (gamma + 1) * (M^2 - 1)))
    b = a - atan(sqrt(M^2 - 1))
    return b

end 

function SolvePrandtlMeyer(nu,gamma)

    F(M)  =  PrandtlMeyer(M,gamma) - nu
    M_out = find_zero(F,2)

    return M_out

end

function TurnAngleLimit(M1,gamma)

    nu_max = pi/2*(((gamma+1)/(gamma-1))^0.5 - 1)
    return nu_max - PrandtlMeyer(M1)

end


function MachAngle(M)
    return asin(1/M)
end

function M2(M1,gamma,TurnAngle)

    nu1 = PrandtlMeyer(M1,gamma)
    nu2 = TurnAngle + nu1
    return  SolvePrandtlMeyer(nu2,gamma)

end


function ApplyExpansionFan(flow,TurnAngle)

    m1    = flow.MachNumber
    gamma = flow.FluidObj.Ratio_of_Specific_Heats

    m2   = M2(m1,gamma,TurnAngle)

    ExpFan_Ratios = Isentropic.Generate_IsentropicRatios(m1,m2,flow.FluidObj.Ratio_of_Specific_Heats)
    t2   = ExpFan_Ratios.Temperature_Ratio  *flow.Temperature
    p2   = ExpFan_Ratios.Pressure_Ratio     *flow.Pressure
    rho2 = ExpFan_Ratios.Density_Ratio      *flow.Density

    ds_fluid = BaseFluid.MakeFluid(flow.FluidObj.Name,flow.FluidObj.Ratio_of_Specific_Heats,flow.FluidObj.Gas_Constant)
    downstream_flow = BaseFluid.MakeFlow(ds_fluid;MachNumber=m2,Temperature=t2,Pressure=p2,Density=rho2)

end


function ApplyExpansionFan_Wrapped(flow,TurnAngle)

    m1    = flow.MachNumber
    gamma = flow.FluidObj.Ratio_of_Specific_Heats

    m2   = M2(m1,gamma,TurnAngle)

    ExpFan_Ratios = Isentropic.Generate_IsentropicRatios(m1,m2,flow.FluidObj.Ratio_of_Specific_Heats)
    t2   = ExpFan_Ratios.Temperature_Ratio  *flow.Temperature
    p2   = ExpFan_Ratios.Pressure_Ratio     *flow.Pressure
    rho2 = ExpFan_Ratios.Density_Ratio      *flow.Density

    ds_fluid = BaseFluid.MakeFluid(flow.FluidObj.Name,flow.FluidObj.Ratio_of_Specific_Heats,flow.FluidObj.Gas_Constant)
    downstream_flow = BaseFluid.MakeFlow(ds_fluid;MachNumber=m2,Temperature=t2,Pressure=p2,Density=rho2)

    incoming_angle = MachAngle(m1)
    outgoing_angle = MachAngle(m2)

    return PM_Wrapper(flow,downstream_flow,incoming_angle,outgoing_angle)


end

end



# println(rad2deg(PMExpansion.PrandtlMeyer(50,1.4)))