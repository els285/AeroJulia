

module BaseFluid

    export Fluid,Flow,SonicVel,Mach2Vel,MakeFlow

    mutable struct Fluid
        # The fluid class defines an arbitrary calorically perfect fluid with thermodynamic material properties
        Name::String
        Ratio_of_Specific_Heats::Float64
        Gas_Constant::Float64
    end

    mutable struct Flow
        # The flow `class` defines the state of a fluid in terms of its thermodynamic variables and flow velocity
        FluidObj::Fluid
        MachNumber::Float64
        Temperature::Float64
        Pressure::Float64
        Density::Float64
    end


    function MakeFlow(fluid,MachNumber;kwargs...)

        kwargs_dict = Dict(kwargs)
        R = fluid.Gas_Constant

        if    haskey(kwargs_dict,:Temperature) && haskey(kwargs_dict,:Pressure)
            T = kwargs_dict[:Temperature]
            p = kwargs_dict[:Pressure]
            rho = p/(R*T)

        elseif haskey(kwargs_dict,:Temperature) && haskey(kwargs_dict,:Density)
            T = kwargs_dict[:Temperature]
            rho = kwargs_dict[:Density]
            p = R*rho*T

        elseif haskey(kwargs_dict,:Pressure) && haskey(kwargs_dict,:Density)
            rho = kwargs_dict[:Density]
            p = kwargs_dict[:Pressure]
            T = p/(R*rho)

        end

        Flow(fluid,MachNumber,T,p,rho)

    end


    function SonicVel(fluid::Fluid,T::Float64)
        
        """
        Returns the acoustic velocity (speed of sound) of a particular in m/s given the temperature
        """    
        return (fluid.Ratio_of_Specific_Heats * fluid.Gas_Constant * T)^0.5

    end

    function Mach2Vel(flow::Flow)

        a = SonicVel(flow.FluidObj,flow.Temperature)
        println(a)
        return a*flow.MachNumber

    end

end



