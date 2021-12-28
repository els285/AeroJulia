include("./BaseFluid.jl")
using .BaseFluid

"""
Taylor-Maccoll: What cone angle gave this shock angle?
1. Guess a shock angle and use the oblique shock relations to guess downstream flow conditions
2. Integrate Taylor-Maccoll equations to determine V_r, V_theta
3. Find theta for V_theta = 0. If this matches the cone angle, found the correct conditions
4. Try agian
"""


module ConicalFlow

    include("./Shocks.jl")
    using .Shocks

    using PyCall
    pushfirst!(pyimport("sys")."path","/home/ethan/Documents/Projects/AeroJulia/GasDynamics/src")
    TMpython = pyimport("justTM")


    function velocity_components(V2,theta_shock,deflec_angle)
        V_r = V2*cos(theta_shock - deflec_angle)
        V_theta = V2*sin(theta_shock - deflec_angle)
        return V_r,V_theta 
    end 


    function generate_conical_flow_from_shock_angle(upstream,theta_shock)

        """
        Generate a conical flow field given an upstream flow and a shock angle
        => Solves the Taylor-MacColl equations numerically to derive the flowfield behind the shockwave, 
            from theta_shock to theta_cone, which is also compute
        """
        # Give theta_shock and upstream flow

        # Compute behind using oblique shock
        downstream = Shocks.ApplyObliqueShock(upstream;beta=theta_shock)
        # Compute the deflection angle of the flow
        deflec_angle = Shocks.ObliqueShock.deflection_angle(upstream.MachNumber,theta_shock,upstream.FluidObj.Ratio_of_Specific_Heats) # Compute here
        #Compute the radial and transverse velocity components of the velocity field
        V_r0 , V_theta0 = velocity_components(downstream.Velocity,theta_shock,deflec_angle)
        # Solve the Taylor-MacColl equation numerically for the flowfield and cone angle
        return TMpython.just_TM(theta_shock,[V_r0,V_theta0])        
    end 
end


f1 = BaseFluid.Fluid("air",1.4,287)

V1 = 5.0
theta_shock_deg = 35.7

flow = BaseFluid.MakeFlow(f1;MachNumber=V1,Density=1.125,Pressure=1e5)


a,b = ConicalFlow.generate_conical_flow_from_shock_angle(flow,deg2rad(theta_shock_deg))




