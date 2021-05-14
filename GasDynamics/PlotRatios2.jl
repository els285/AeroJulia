module GasDynamic_Plots

export Unpack_Arguments,Plot_Ratios,PlotRepeat
using Plots

function Unpack_Arguments(fluid;kwargs)

    gamma = fluid.Ratio_of_Specific_Heats
    kwargs = Dict(kwargs)

    """
    Imports custom Mach number range, or uses default
    """
    if haskey(kwargs,:M_range)
        M1_array = [kwargs[:M_range][1]:0.1:kwargs[:M_range][2];]
    else
        M1_array = [1.0:0.1:15;]
    end

    return gamma,M1_array

end

function Plot_TBM_fromDict(D)

    gr()
    plt = plot()
    for (M,v) in D
        plot!(v[2],v[1],label=M,lw=2,legend=:bottomright)
    end

    display(plt)
    readline()

end 


function Plot_Ratios(xdata,ydata_dic)

    xlabel          = keys(xdata)
    ylabels_column  = collect(keys(ydata_dic))
    ylabels         = reshape(ylabels_column,(1,length(ylabels_column)))
    x_axis_data     = collect(values(xdata))
    data            = collect(values(ydata_dic) )

    gr()
    plt = plot(x_axis_data, data,label=ylabels,lw=2,legend=:topleft)


    display(plt)
    readline()

end


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