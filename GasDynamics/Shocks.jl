module Shocks

"""
The Shocks module contains two sub-modules: NormalShock and ObliqueShock
    These sub-modules contain the equations for calculating thermodynamic ratios across the relevant shock wave

The ShockRatio structure stores the calculated thermodynamic ratios

The ApplyShock functions take an input Flow structure,
    and apply the relevant shock relations to generate a new downstream Flow structure
"""

include("./BaseFluid.jl")
using .BaseFluid


export NormalRatios,ApplyNormalShock,ObliqueRatios,ApplyObliqueShock

mutable struct ShockRatio
    """
    Struct for storing shock wave ratios
    """
    Mach_Ratio
    Temperature_Ratio
    Pressure_Ratio
    Density_Ratio
end


module NormalShock

    """
    Normal Shock Module - Gas Dynamics

    Ethan Simpson - 8th April 2021

    Selection of functions for generating thermodynamic and MachNumber parameters downstream of some normal shock,
        for some specified upwind flow condition
    """
        
    function M2(M1,gamma)

        return (  (M1^2*(gamma-1) + 2)  /  (2*gamma*M1^2 - (gamma-1)) )^0.5

    end

    function P2(M1,gamma)

        return (2*gamma*M1^2 - (gamma -1)) / (gamma+1)

    end

    function T2(M1,gamma)

        return (2*gamma*M1^2 - (gamma-1))*((gamma-1)*M1^2+2) / ((gamma+1)^2*M1^2)

    end

    function Rho2(M1,gamma)

        return ((gamma+1)*M1^2)/((gamma-1)*M1^2+2)

    end 

end

function NormalRatios(MachNumber,gamma)

    Mach_Ratio          = NormalShock.M2(MachNumber,gamma)/MachNumber
    Temperature_Ratio   = NormalShock.T2(MachNumber,gamma)
    Pressure_Ratio      = NormalShock.P2(MachNumber,gamma)
    Density_Ratio       = NormalShock.Rho2(MachNumber,gamma)

    return ShockRatio(Mach_Ratio,Temperature_Ratio,Pressure_Ratio,Density_Ratio)

end

function ApplyNormalShock(flow)

    NormalShockRatio = NormalRatios(flow.MachNumber,flow.FluidObj.Ratio_of_Specific_Heats)
    m2   = NormalShockRatio.Mach_Ratio         *flow.MachNumber
    t2   = NormalShockRatio.Temperature_Ratio  *flow.Temperature
    p2   = NormalShockRatio.Pressure_Ratio     *flow.Pressure
    rho2 = NormalShockRatio.Density_Ratio      *flow.Density

    ds_fluid = BaseFluid.Fluid(flow.FluidObj.Name,flow.FluidObj.Ratio_of_Specific_Heats,flow.FluidObj.Gas_Constant)
    downstream_flow = BaseFluid.Flow(ds_fluid,m2,t2,p2,rho2)

    return downstream_flow

end


##############################################################################################################

module ObliqueShock

    """
    Oblique Shock Module - Gas Dynamics

    Ethan Simpson - 8th April 2021

    Selection of functions for generating thermodynamic and MachNumber parameters downstream of some oblique shock,
        for some specified upwind flow condition
    """

    function ConeAngle(M1,beta,gamma)

        # Gives cone angle for input Mach number and shock angle
        tan_theta = 2*cot(beta)*(M1^2*sin(beta)^2-1)/(M1^2*(gamma + cos(2*beta))+2)

        return atan(tan_theta)

    end


    function M2(M1,beta,gamma)

        theta = ConeAngle(M1,beta,gamma)

        return 1/(sin(beta - theta))*((1 + 0.5*(gamma-1)*M1^2*sin(beta)^2)/(gamma*M1^2*sin(beta)^2) - 0.5*(gamma-1))^0.5

    end


    function P2(M1,beta,gamma)

        return 1 + 2*gamma/(gamma+1) *(M1^2*sin(beta)^2 - 1)

    end

    function T2(M1,beta,gamma)

        return (2*gamma*M1^2*sin(beta)^2 - (gamma-1)) * ((gamma-1)*M1^2*sin(beta)^2 +2) / ((gamma+1)^2*M1^2*sin(beta)^2)

    end

    function Rho2(M1,beta,gamma)

        return ((gamma+1)*M1^2*sin(beta)^2)  / ((gamma-1)*M1^2*sin(beta)^2 + 2 )

    end 


end

    ##########################

function ObliqueRatios(MachNumber,gamma,beta)

    Mach_Ratio          = ObliqueShock.M2(MachNumber,beta,gamma)/MachNumber
    Temperature_Ratio   = ObliqueShock.T2(MachNumber,beta,gamma)
    Pressure_Ratio      = ObliqueShock.P2(MachNumber,beta,gamma)
    Density_Ratio       = ObliqueShock.Rho2(MachNumber,beta,gamma)

    return ShockRatio(Mach_Ratio,Temperature_Ratio,Pressure_Ratio,Density_Ratio)

end

function ApplyObliqueShock(flow,beta)

    ObliqueShockRatio = ObliqueRatios(flow.MachNumber,flow.FluidObj.Ratio_of_Specific_Heats,beta)
    m2   = ObliqueShockRatio.Mach_Ratio         *flow.MachNumber
    t2   = ObliqueShockRatio.Temperature_Ratio  *flow.Temperature
    p2   = ObliqueShockRatio.Pressure_Ratio     *flow.Pressure
    rho2 = ObliqueShockRatio.Density_Ratio      *flow.Density
    
    ds_fluid = BaseFluid.Fluid(flow.FluidObj.Name,flow.FluidObj.Ratio_of_Specific_Heats,flow.FluidObj.Gas_Constant)
    downstream_flow = BaseFluid.Flow(ds_fluid,m2,t2,p2,rho2)

    return downstream_flow

end


end