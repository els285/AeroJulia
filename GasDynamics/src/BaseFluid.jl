
module BaseFluid

    export Fluid,Flow,SonicVel,Mach2Vel,MakeFlow,MakeFluid,Air

    mutable struct Fluid
        # The fluid class defines an arbitrary calorically perfect fluid with thermodynamic material properties
        Name::String
        Ratio_of_Specific_Heats::Float64
        Gas_Constant::Float64
    end

    function MakeFluid(Name,Ratio_of_Specific_Heats,Gas_Constant)
        return Fluid(Name,Ratio_of_Specific_Heats,Gas_Constant)
    end

    mutable struct Flow
        # The flow `class` defines the state of a fluid in terms of its thermodynamic variables and flow velocity
        FluidObj::Fluid
        MachNumber::Float64
        Velocity::Float64
        Temperature::Float64
        Pressure::Float64
        Density::Float64
    end

    function SonicVel(fluid::Fluid,T::Float64)
        
        """
        Returns the acoustic velocity (speed of sound) of a particular in m/s given the temperature
        """    
        return (fluid.Ratio_of_Specific_Heats * fluid.Gas_Constant * T)^0.5

    end

    function Mach2Vel(MachNumber,Temperature,flu_obj)

        """ Taking MachNumber, temperature and FluidObj, computes acoustic velocity and returns velocity """

        a = SonicVel(flu_obj,Temperature)
        return a*MachNumber

    end


    function Vel2Mach(Velocity,Temperature,flu_obj)

        """ Taking Velocity, temerpature and FluidObj, computes acoustic velocity and returns MachNumber """

        a = SonicVel(flu_obj,Temperature)
        return Velocity/a

    end


    function MakeFlow(fluid;kwargs...)

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


        if     haskey(kwargs_dict,:MachNumber) && !haskey(kwargs_dict,:Velocity)
            MachNumber = kwargs_dict[:MachNumber]
            Velocity   = Mach2Vel(MachNumber,T,fluid)
        elseif haskey(kwargs_dict,:Velocity)   && !haskey(kwargs_dict,:MachNumber)
            Velocity   = kwargs_dict[:Velocity]
            MachNumber = Vel2Mach(Velocity,T,fluid)
        end

        Flow(fluid,MachNumber,Velocity,T,p,rho)

    end

    # Common fluids

    Air = MakeFluid("Standard dry air @ 293K",1.4,287) # Standard dry air 293K, useful for most problems

end



