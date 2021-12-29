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
    profile_points
    internal_angles
end


function gen_body(profile_points)
    tup_profile = []
    for p in profile_points
        push!(tup_profile,(x=p[1],y=p[2]))
    end 
    internal_angles = generate_internal_angles(tup_profile)
    return Body(tup_profile,internal_angles)
end




# function generate_body(profile_points)
#     reflected_points = []#Array(Int,size(profile_points))
#     for point in profile_points
#         push!(reflected_points,(point[1],-point[2]))
#     end 
#     profile = Set(vcat(profile_points,reflected_points))
#     return Body(profile)
# end 



function get_tan_angle(p1,p2)
    dy = (p2.y - p1.y)
    dx = (p2.x - p1.x)
    return atan(dy,dx)
end


function generate_internal_angles(profile)

    internal_angles = []
    a1 = get_tan_angle(profile[1],profile[2])
    push!(internal_angles,a1)

    N = length(profile)

    for i in 1:(N-2)
        p1 = profile[i]
        p2 = profile[i+1]
        p3 = profile[i+2]
        t1 = get_tan_angle(p1,p2)
        t2 = get_tan_angle(p2,p3)
        theta = t2 -t1
        push!(internal_angles,theta)
    end 

    aN = get_tan_angle(profile[N-1],profile[N])
    push!(internal_angles,aN)

    return internal_angles

end


function apply(Upstream_flow, Body)

#     # Oblique shocks and PM expansions
    list_of_flow_conditions = []
    push!(list_of_flow_conditions,Upstream_flow)
    current_flow = Upstream_flow

    for angle in Body.internal_angles[1:length(Body.internal_angles)-1]
        if angle > 0
            new_flow = ApplyObliqueShock(current_flow;beta=angle)
            # Apply compression oblique shock
        elseif angle < 0 
            new_flow = ApplyExpansionFan(current_flow,angle)
            # Apply PM expansion
        else 
            new_flow = current_flow
        end 
        push!(list_of_flow_conditions,new_flow)
    end

    # The final flow change is a trailing-edge shock wave of same angle
    fin_flow = ApplyObliqueShock(last(list_of_flow_conditions),beta=Body.internal_angles[1])
    push!(list_of_flow_conditions,fin_flow)
    return list_of_flow_conditions

end
        




prof_points = [(0,0),(10,3),(20,0)]
b1 = gen_body(prof_points)
f1 = BaseFluid.MakeFlow(Air;MachNumber=10,Density=1.125,Pressure=1e5)

flow_solution = apply(f1,b1 )

for f in flow_solution

    println(f)

end

