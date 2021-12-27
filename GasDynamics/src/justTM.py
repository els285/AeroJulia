import numpy as np 

from scipy.interpolate import interp1d,UnivariateSpline
from scipy.integrate import odeint
from mpmath import cot

pi =np.pi


def TM_fromMATLAB(g,t):
    gamma=1.4
    a = (gamma-1)/2
    Numerator = a*(2*((g[0])**3-g[0])+(g[1]-g[1]*g[0]**2-g[1]**3)*cot(t)-2*g[0]*g[1]**2)-(g[0]**2*g[1]**2)
    Denominator = a*((g[1])**2+(g[0])**2 - 1) + (g[1])**2
    return [g[1],Numerator/Denominator]
    # du[0] = g[1]
    # du[1] = Numerator/Denominator


def taylor_maccoll(u,t):

    x , y = u

    gamma = 1.4
    P = (gamma-1)*0.5*(1 - x**2 - y**2)
    du = [y,(-P*(2*x + y*cot(t)) +x*y**2)/(P - y**2)]
    return du


# def taylor_maccoll(y, theta, gamma=1.4):
#     # Taylor-Maccoll function
#     # Source: https://www.grc.nasa.gov/www/k-12/airplane/coneflow.html
#     v_r, v_theta = y
#     dydt = [
#         v_theta,
#         (v_theta ** 2 * v_r - (gamma - 1) / 2 * (1 - v_r ** 2 - v_theta ** 2) * (2 * v_r + v_theta / np.tan(theta))) / ((gamma - 1) / 2 * (1 - v_r ** 2 - v_theta ** 2) - v_theta ** 2) 
#     ]
#     return dydt


def just_TM(beta_angle,y0):

    # Theta range
    thetas = np.linspace(-beta_angle, -5*pi/180, 500)
    thetas_array = thetas.reshape((len(thetas),1))

    # Define the ODE solution
    sol = odeint(taylor_maccoll, y0, thetas)


    # Define the flowfield
    flow_field = np.concatenate((thetas_array,sol),axis=1)

    idx = (np.abs(flow_field[:,2] - 0.0)).argmin()
    physical_flow_field = flow_field[:idx,:]

    # Compute the 
    freduced = UnivariateSpline(flow_field[:,0], flow_field[:,2], s=0)
    cone_angle = -freduced.roots()[0]
    
    return physical_flow_field, cone_angle






