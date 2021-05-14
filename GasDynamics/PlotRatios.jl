module GasDynamic_Plots

export Make_Ratio_Plot
using Plots

include("Shocks.jl")
using .Shocks

include("Isentropic.jl")
using .Isentropic


function Unpack_Arguments(fluid;kwargs)

    gamma = fluid.Ratio_of_Specific_Heats

    """
    Imports custom Mach number range, or uses default
    """
    if haskey(kwargs,:M_range)
        M1_array = [kwargs[:M_range][1]:0.1:kwargs[:M_range][2];]
    else
        M1_array = [1.0:0.1:8;]
    end

    return gamma,M1_array

end


function Make_Ratio_Plot(fluid,pheno;kwargs...)



    """Fluid Properties"""
    gamma = fluid.Ratio_of_Specific_Heats



    """
    Determines type of gas dynamic phenomenon to plot
    """
    kwargs = Dict(kwargs)
    if pheno    == "Normal"
        Ratios_Object = Shocks.NormalRatios
        input = [gamma]

    elseif pheno == "Oblique"
        Ratios_Object = Shocks.ObliqueRatios
        @assert haskey(kwargs,:beta) "Beta angle must be provided"
        # print(kwargs(:beta))
        # readlline()
        input = [gamma,kwargs[:beta]]

    elseif pheno == "Isentropic"
        Ratios_Object = Isentropic.Generate_IsentropicRatios
        @assert haskey(kwargs,:M2) "M2 (downstream Mach Number) must be specified"
        input = [kwargs[:M2],gamma]

    end

    """
    Imports custom Mach number range, or uses default
    """
    if haskey(kwargs,:M_range)
        M1_array = [kwargs[:M_range][1]:0.1:kwargs[:M_range][2];]
    else
        M1_array = [1.0:0.1:8;]
    end

    # print(input)
    # readline()
    M2_array     = Float64[]
    Temperatures = Float64[]
    Pressures    = Float64[]
    Densities    = Float64[]

    for M1 in M1_array
        shock_conditions = Ratios_Object(M1,input...)
        # push!(M2_array,shock_conditions.Mach_Ratio)
        push!(Temperatures,shock_conditions.Temperature_Ratio)
        push!(Pressures,shock_conditions.Pressure_Ratio)
        push!(Densities,shock_conditions.Density_Ratio)
    end


    """
    Plot
    """
    gr()
    plt = plot(M1_array, [Temperatures,Pressures,Densities],
                label=["Mach Ratio" "Temperature" "Pressure" "Density"],    
                lw=2)

    display(plt)
    readline()



end



# function Plot_Ratios(data)

#     M1_array     = data.M1
#     M2_array     = data.M_ratios
#     Temperatures = data.T_ratios
#     Pressures    = data.P_ratios
#     Densities    = data.Rho_ratios

#     gr()
#     plt = plot(M1_array, [M_ratio_array,Temperatures,Pressures,Densities],label=["Mach Ratio" "Temperature" "Pressure" "Density"],lw=2)

#     display(plt)
#     readline()

# end


# function Plot_All_Ratios(Fluid,Ratios_Object;kwargs...)

#     # if string(Ratios_Object)
#     print(Ratios_Object isa Main.GasDynamics.Shocks.NormalRatios)
#     # print(typeof(Ratios_Object))

#     # print(string(Ratios_Object))
#     readline()

    

#     gr()
#     gamma = Fluid.Ratio_of_Specific_Heats

#     M1_array = [1.0:0.1:10;]

#     M2_array     = Float64[]
#     Temperatures = Float64[]
#     Pressures    = Float64[]
#     Densities    = Float64[]

#     for M1 in M1_array
#         shock_conditions = Ratios_Object(M1,gamma)
#         push!(M2_array,shock_conditions.Mach_Ratio)
#         push!(Temperatures,shock_conditions.Temperature_Ratio)
#         push!(Pressures,shock_conditions.Pressure_Ratio)
#         push!(Densities,shock_conditions.Density_Ratio)
#     end


#     plt = plot(M1_array, [M2_array,Temperatures,Pressures,Densities],label=["Mach Ratio" "Temperature" "Pressure" "Density"],lw=2)

#     display(plt)
#     readline()

# end


end