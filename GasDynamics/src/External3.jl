# External flow
include("./Shocks.jl")
using .Shocks


include("./BaseFluid.jl")
using .BaseFluid

include("./PMExpansion.jl")
using .PMExpansion

"""
Script for external flow simulation of horizontally-symmetric bodies at zero angle of attack in supersonic flows
Quasi-one-dimensional flow
"""


mutable struct Body
    profiles
    angles
    aoa
end


function gen_body(upper_profile_points,lower_profile_points,aoa)
    
    # Re-write as NamedTuple coordinates
    upper_tup_profile = [(x=p[1],y=p[2]) for p in upper_profile_points]
    lower_tup_profile = [(x=p[1],y=p[2]) for p in lower_profile_points]
    
    # Compute angles through which flow rotates
    upper_internal_angles = generate_flowchange_angles(upper_tup_profile,aoa)#,"upper")
    lower_internal_angles = generate_flowchange_angles(lower_tup_profile,aoa)#,"lower")
    
    # Assign to NamedTuples
    profiles = (upper=upper_tup_profile,lower=lower_tup_profile)
    angles   = (upper=upper_internal_angles,lower=lower_internal_angles)
    return Body(profiles,angles,aoa)
end


function compute_aero_forces(Body,flow)
    # pass
end



function get_tan_angle(p1,p2)
    dy = (p2.y - p1.y)
    dx = (p2.x - p1.x)
    return atan(dy,dx)
end


function generate_flowchange_angles(profile,aoa)

    internal_angles = []
    N = length(profile)
    
    a1 =   aoa - get_tan_angle(profile[1],profile[2])
    aN =  get_tan_angle(profile[N-1],profile[N]) - aoa

    # if edge=="upper"
    #     a1 = aoa - get_tan_angle(profile[1],profile[2])
    #     aN = get_tan_angle(profile[N-1],profile[N]) - aoa
    # elseif edge=="lower"
    #     a1 =   aoa - get_tan_angle(profile[1],profile[2])
    #     aN =  get_tan_angle(profile[N-1],profile[N]) - aoa
    # end
        
    push!(internal_angles,a1)

    for i in 1:(N-2)
        p1 = profile[i]
        p2 = profile[i+1]
        p3 = profile[i+2]
        t1 = get_tan_angle(p1,p2)
        t2 = get_tan_angle(p2,p3)
        theta = t1-t2
        push!(internal_angles,theta)
    end 

    push!(internal_angles,aN)

    return internal_angles

end


function apply(Upstream_flow, Body)

    """
    Takes Upstream flow which defines the incident flow velocity and theromdynamic state
    Takes Body which contains the geometry of the upper and lower surfaces

    Flow assumed to be left -> right
    The definition of angles is such that:
        For Upper-Surface, t1 - t2 < 0 implies that the flow moves "left" through compression shock
                           t1 - t2 > 0 implies that the flow moves "right" through expansion fan
        For Lower-Surace, this is reversed.
    """

    # UPPER EDGE
    upper_flows = []
    push!(upper_flows,Upstream_flow)
    current_flow = Upstream_flow

    U = Body.angles.upper
    for angle in U[1:length(U)]
        if angle < 0
            # Apply compression oblique shock
            output = ApplyObliqueShock_Wrapped(current_flow;theta=abs(angle))
            new_flow = output.Downstream
        elseif angle > 0 
            # Apply PM expansion
            output = ApplyExpansionFan_Wrapped(current_flow,angle)
            new_flow = output.Downstream
        else 
            output = nothing
            new_flow = current_flow
        end 
        current_flow = new_flow
        push!(upper_flows,output)
    end
   
    # The final flow change is a trailing-edge shock wave of same angle
    # fin_flow = ApplyObliqueShock(last(upper_flows),beta=Body.internal_angles[1])
    # push!(upper_flows,fin_flow)  
    

    ## DOING LOWER EDGE

    lower_flows = []
    push!(lower_flows,Upstream_flow)
    current_flow = Upstream_flow    

    L = Body.angles.lower
    for angle in L[1:length(L)]
        if angle > 0
            # Apply compression oblique shock
            output = ApplyObliqueShock_Wrapped(current_flow;theta=angle)
            new_flow = output.Downstream

        elseif angle < 0 
            # Apply PM expansion
            output = ApplyExpansionFan_Wrapped(current_flow,abs(angle))
            new_flow = output.Downstream
        else 
            output = nothing
            new_flow = current_flow
        end 
        current_flow = new_flow
        push!(lower_flows,output)
    end

    # fin_flow = ApplyObliqueShock(last(lower_flows),beta=lastBody.internal_angles)
    # push!(lower_flows,fin_flow)  

    return (upper = upper_flows , lower = lower_flows)

end


# function plot_shock_expansion


        


upper_points = [(0,0),(10,1),(20,0)]
lower_points = [(0,0),(10,-1),(20,0)]
aoa = deg2rad(0)
b1 = gen_body(upper_points,lower_points,aoa)


f1 = BaseFluid.MakeFlow(Air;MachNumber=5,Density=1.125,Pressure=1e5)

flow_solution = apply(f1,b1)


using JSON

open("foo.json","w") do f
    d = Dict("Body" => b1 , "FlowSolution" => flow_solution)
    JSON.print(f,d)
    # JSON.print(f, b1)
    # JSON.print(f, flow_solution)
end


# for f in flow_solution.upper
#     println(f.MachNumber)
# end


# for f in flow_solution.lower
#     println(f.MachNumber)
# end