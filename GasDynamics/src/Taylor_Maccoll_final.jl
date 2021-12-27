# __precompile__()
using DifferentialEquations
# using OrdinaryDiffEq
# Taylor-Maccoll Analysis

include("./Shocks.jl")
include("./BaseFluid.jl")
using .Shocks
using .BaseFluid

"""
Taylor-Maccoll: What cone angle gave this shock angle?
1. Guess a shock angle and use the oblique shock relations to guess downstream flow conditions
2. Integrate Taylor-Maccoll equations to determine V_r, V_theta
3. Find theta for V_theta = 0. If this matches the cone angle, found the correct conditions
4. Try agian
"""

function TM_fromMATLAB(du,g,p,t)
    gamma=1.4
    a = (gamma-1)/2
    Numerator = a*(2*((g[1])^3-g[1])+(g[2]-g[2]*g[1]^2-g[2]^3)*cot(t)-2*g[1]*g[2]^2)-(g[1]^2*g[2]^2)
    Denominator = a*((g[2])^2+(g[1])^2 - 1) + (g[2])^2
    du[1] = g[2]
    du[2] = Numerator/Denominator
end

# % Numerator and denominator for zdot calculation below

# function f(du,z,p,theta)
#     gamma=1.4
#     A = (gamma-1)/2
#     num = (-2*A*z[1]) - (A*z[2]*cot(theta)) + (2*A*z[1]^3) + (A*z[1]^2*z[2]*cot(theta)) + (2*A*z[1]*z[2]^2) + (A*z[2]^3*cot(theta)) + (z[1]*z[2]^2)
#     den = A*(1-z[1]^2-z[2]^2) - z[2]^2
#     du[1] = z[2]
#     du[2] = num/den
# end

# function f(du,u,p,t)

#     x , y = u

#     du[1] = y    # dv_r/dtheta = V_theta 

#     gamma = 1.4
#     P = (gamma-1)*0.5*(1 - x^2 - y^2)
#     du[2] = (-P*(2*x + y*cot(t)) +x*y^2)/(P - y^2)
#     # du[2] = x

# end

# function P(u,p,t) 

#     V = u
#     theta = t

#     output = zeros(1, 2)

#     output[1] = V[2]

#     gamma = 1.4

#     P = (gamma-1)*0.5*(1 - V[1]^2 - V[2]^2)

#     output[2] = (-P*(2*V[1] + V[2]*cot(theta)) +V[1]*V[2]^2)/(P - V[2]^2)

#     return output

# end

# def taylor_maccoll(y, theta, gamma=1.4):
#     # Taylor-Maccoll function
#     # Source: https://www.grc.nasa.gov/www/k-12/airplane/coneflow.html
#     v_r, v_theta = y
#     dydt = [
#         v_theta,
#         (v_theta ^ 2 * v_r - (gamma - 1) / 2 * (1 - v_r ^ 2 - v_theta ** 2) * (2 * v_r + v_theta / np.tan(theta))) / ((gamma - 1) / 2 * (1 - v_r ** 2 - v_theta ** 2) - v_theta ** 2) 
#     ]
#     return dydt


# function f(du,u,p,theta)

#     v_r,v_theta = u

#     du[1] = v_theta

#     du[2] = (v_theta ^ 2 * v_r - (gamma - 1) / 2 * (1 - v_r ^ 2 - v_theta ^ 2) * (2 * v_r + v_theta / tan(theta))) / ((gamma - 1) / 2 * (1 - v_r ^ 2 - v_theta ^ 2) - v_theta ^ 2) 

#     return du 

# end

function deflection_angle(M_free,theta_shock)
    gamma = 1.4
    tandfa = 2/tan(theta_shock)*(M_free^2*sin(theta_shock)^2-1)/(M_free^2*(gamma+cos(2*theta_shock))+2)
    return atan(tandfa)
end


# println(deflection_angle(2,0.5))

gamma =1.4

function non_dimensionalised_velocity(M2)

    return (2/((gamma-1)*M2^2)+1)^(-0.5)

end



f1 = BaseFluid.Fluid("air",1.4,287)

flow = BaseFluid.MakeFlow(f1;MachNumber=3.0,Density=1.125,Pressure=1e5)


# M1 = 2.1
theta_shock = 31.7*pi/180
ds = Shocks.ApplyObliqueShock(flow;beta=0.5)

V_prime = non_dimensionalised_velocity(ds.MachNumber)


deflec_angle = deflection_angle(flow.MachNumber,theta_shock)
V2 = ds.Velocity

V_r0 = V_prime*cos(theta_shock - deflec_angle)
V_theta0 = V_prime*sin(theta_shock - deflec_angle)

println(V_r0)
println(V_theta0)

theta_span = (theta_shock,0.1*pi/180)

prob = ODEProblem(TM_fromMATLAB,[V_r0,V_theta0],theta_span)

num_sol = solve(prob,Tsit5(),reltol=1e-6,abstol=1e-6)

println(last(num_sol[1]))
println(last(num_sol[2]))

