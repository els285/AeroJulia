module Shocks

using Roots


"""
The Shocks module contains two sub-modules: NormalShock and ObliqueShock
    These sub-modules contain the equations for calculating thermodynamic ratios across the relevant shock wave

The ShockRatio structure stores the calculated thermodynamic ratios

The ApplyShock functions take an input Flow structure,
    and apply the relevant shock relations to generate a new downstream Flow structure

The plotting functions slow down the whole code because the Plots package is slow
Think about upgrading to a new plotting package
"""

include("./BaseFluid.jl")
using .BaseFluid

include("GDPlots.jl")
using .GasDynamic_Plots



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

function ApplyNormalShock(input_flow)

    """
    ApplyNormalShock computes a ShockRatio object for a given input_flow 
    """

    NormalShockRatio = NormalRatios(input_flow.MachNumber,input_flow.FluidObj.Ratio_of_Specific_Heats)
    m2   = NormalShockRatio.Mach_Ratio         *input_flow.MachNumber
    t2   = NormalShockRatio.Temperature_Ratio  *input_flow.Temperature
    p2   = NormalShockRatio.Pressure_Ratio     *input_flow.Pressure
    rho2 = NormalShockRatio.Density_Ratio      *input_flow.Density

    # ds_fluid = input_flow.FluidObj#BaseFluid.Fluid(flow.FluidObj.Name,flow.FluidObj.Ratio_of_Specific_Heats,flow.FluidObj.Gas_Constant) 
    downstream_flow = BaseFluid.Flow(input_flow.FluidObj,m2,t2,p2,rho2)

    return downstream_flow

end


##############################################################################################################

module ObliqueShock

    using Roots

    """
    Oblique Shock Module - Gas Dynamics

    Ethan Simpson - 8th April 2021

    Selection of functions for generating thermodynamic and MachNumber parameters downstream of some oblique shock,
        for some specified upwind flow condition


    Oblique shocks can lead to supersonic (weak shock) or subsonic (strong shock) downstream conditions
    """

    function ConeAngle(M1,beta,gamma)
        # Gives cone angle for input Mach number and shock angle
        tan_theta = 2*cot(beta)*(M1^2*sin(beta)^2-1)/(M1^2*(gamma + cos(2*beta))+2)
        return atan(tan_theta)
        end


    function ShockAngle(M1,theta,gamma)

        # Cone angle returns the angle theta
        # Use Newton-Raphson to solve for theta_out
        # Should be two solutions; this method returns the weak solution

        F(beta) = ConeAngle(M1,beta,gamma) - theta
        beta_out = Roots.find_zero(F,1)
        return beta_out

    end 


    function Minimum_ShockAngle(M1,gamma)
        # The minimum shock angle is limited by Mach number: fow low Mach numbers, even zero cone angle will produce a non-zero beta               
        return ShockAngle(M1,0,gamma)
    end




    function M2(M1,beta,gamma)

        theta = ConeAngle(M1,beta,gamma)

        numer = 1 + 0.5*(gamma-1)*M1^2*(sin(beta))^2
        denom = gamma*M1^2*sin(beta)^2 - 0.5*(gamma-1)

        return 1/(sin(beta - theta))*(numer/denom)^0.5

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


    function deflection_angle(M_free,theta_shock,gamma)
        tandfa = 2/tan(theta_shock)*(M_free^2*sin(theta_shock)^2-1)/(M_free^2*(gamma+cos(2*theta_shock))+2)
        return atan(tandfa)
    end



end



function non_dimensionalised_velocity(M2,gamma)
    return (2/((gamma-1)*M2^2)+1)^(-0.5)
end

    ##########################

function ObliqueRatios(MachNumber,gamma,beta)

    Mach_Ratio          = ObliqueShock.M2(MachNumber,beta,gamma)/MachNumber
    Temperature_Ratio   = ObliqueShock.T2(MachNumber,beta,gamma)
    Pressure_Ratio      = ObliqueShock.P2(MachNumber,beta,gamma)
    Density_Ratio       = ObliqueShock.Rho2(MachNumber,beta,gamma)

    return ShockRatio(Mach_Ratio,Temperature_Ratio,Pressure_Ratio,Density_Ratio)

end

function ApplyObliqueShock(flow;angles...)

    # An ObliqueShock can be initialised with either theta or beta
    if      haskey(angles,:beta) && !(haskey(angles,:theta))
        beta  = angles[:beta]
        theta = ObliqueShock.ConeAngle(flow.MachNumber,beta,flow.FluidObj.Ratio_of_Specific_Heats)
    elseif  haskey(angles,:theta) && !(haskey(angles,:beta))
        theta = angles[:theta]
        beta  = ObliqueShock.ShockAngle(flow.MachNumber,theta,flow.FluidObj.Ratio_of_Specific_Heats)
    end

    # Generate the oblique shock ratios
    ObliqueShockRatio = ObliqueRatios(flow.MachNumber,flow.FluidObj.Ratio_of_Specific_Heats,beta)

    # Apply to original flow
    m2   = ObliqueShockRatio.Mach_Ratio         *flow.MachNumber
    t2   = ObliqueShockRatio.Temperature_Ratio  *flow.Temperature
    p2   = ObliqueShockRatio.Pressure_Ratio     *flow.Pressure
    rho2 = ObliqueShockRatio.Density_Ratio      *flow.Density
    
    # Construct new flow - why do we need to make a new fluid?
    ds_fluid = BaseFluid.MakeFluid(flow.FluidObj.Name,flow.FluidObj.Ratio_of_Specific_Heats,flow.FluidObj.Gas_Constant)
    downstream_flow = BaseFluid.MakeFlow(ds_fluid;MachNumber=m2,Temperature=t2,Pressure=p2)

    return downstream_flow

end

function Plot_NormalRatios(fluid;kwargs...)

    gamma,M1_array = Unpack_Arguments(fluid;kwargs)
    M1_array = [1.0:0.1:10;]

    M2_array     = Float64[]
    Temperatures = Float64[]
    Pressures    = Float64[]
    Densities    = Float64[]

    for M1 in M1_array
        shock_conditions = NormalRatios(M1,gamma)
        push!(M2_array,shock_conditions.Mach_Ratio)
        push!(Temperatures,shock_conditions.Temperature_Ratio)
        push!(Pressures,shock_conditions.Pressure_Ratio)
        push!(Densities,shock_conditions.Density_Ratio)
    end

    data = Dict([("M2/M1",M2_array),
                 ("T2/T1",Temperatures),
                 ("P2/P1",Pressures),
                 ("rho1/rho2 ",Densities)])

    Plot_Ratios(Dict([("M1",M1_array)]) , data)

end


function Plot_ObliqueRatios_M1(fluid,beta;kwargs...)

    """
    Plot shows shock ratios as a function of upstream Mach number
    """

    gamma,M1_array = Unpack_Arguments(fluid;kwargs)

    M2_array     = Float64[]
    Temperatures = Float64[]
    Pressures    = Float64[]
    Densities    = Float64[]

    for M1 in M1_array
        shock_conditions = ObliqueRatios(M1,gamma,beta)
        push!(M2_array,shock_conditions.Mach_Ratio)
        push!(Temperatures,shock_conditions.Temperature_Ratio)
        push!(Pressures,shock_conditions.Pressure_Ratio)
        push!(Densities,shock_conditions.Density_Ratio)
    end

    data = Dict([("M2/M1",M2_array),
    ("T2/T1",Temperatures),
    ("P2/P1",Pressures),
    ("rho1/rho2 ",Densities)])

    Plot_Ratios(Dict([("M1",M1_array)]) , data)


end


function Plot_ObliqueRatios_beta(fluid,M1;kwargs...)

    """
    Show shock ratios as a function of shock angle
    """

    gamma = fluid.Ratio_of_Specific_Heats

    beta_min = ObliqueShock.Minimum_ShockAngle(M1,gamma)

    beta_array = [beta_min:0.01:pi/2;]


    M2_array     = Float64[]
    Temperatures = Float64[]
    Pressures    = Float64[]
    Densities    = Float64[]

    for beta in beta_array
        shock_conditions = ObliqueRatios(M1,gamma,beta)
        push!(M2_array,shock_conditions.Mach_Ratio)
        push!(Temperatures,shock_conditions.Temperature_Ratio)
        push!(Pressures,shock_conditions.Pressure_Ratio)
        push!(Densities,shock_conditions.Density_Ratio)
    end

    data = Dict([("M2/M1",M2_array),
    ("T2/T1",Temperatures),
    ("P2/P1",Pressures),
    ("rho1/rho2 ",Densities)])

    Plot_Ratios(Dict([("beta",beta_array*180/pi)]) , data)

end



function Plot_ThetaBetaM(fluid;kwargs...)

    """
    Generating the famous theta-beta-M plot
    """

    gamma = fluid.Ratio_of_Specific_Heats

    if haskey(kwargs,:M_range)
        M_range = kwargs(:M_range)
    else
        M_range = [1.1,1.2,1.4,1.6,1.8,2.0,2.25,2.5,3.0,4.0,5.0,8.0,10.0,15.0]
    end

    D = Dict()


    for M in M_range

        beta_min    = ObliqueShock.Minimum_ShockAngle(M,gamma)
        beta_array  = collect(LinRange(beta_min,pi/2,100) )
        theta_array = Float64[]

        for beta in beta_array
            push!(theta_array,ObliqueShock.ConeAngle(M,beta,gamma))
        end 

        D[M] = (beta_array,theta_array)

    end

    Plot_TBM_fromDict(D)


end



end