# Fanno flow
"""
Fanno flow is irreversible, adiabatic flow with friction
Fanno flow occurs in a duct when area and mass flow are held constant
For steady quasi-one-dimensional flows
"""


module Fanno

    function f(M,gamma)
        return (2/(gamma+1)*(1+(gamma-1)/2*M^2))
    end

    function P_ratio(M,gamma)
        return 1/M *f(M,gamma)^(-0.5)
    end 

    function Rho_ratio(M,gamma)
        return 1/M*f(M,gamma)^(0.5)
    end

    function T_ratio(M,gamma)
        return f(M,gamma)^(-1)
    end 

    function V_ratio(M,gamma)
        return M*f(M,gamma)^(-0.5)
    end

    function P_stag_ratio(M,gamma)
        return 1/M*f(M,gamma)^((gamma+1)/(2*(gamma-1)))
    end 

    function Fanno_parameter(M,gamma)
        return (1-M^2)/(gamma*M^2) + (gamma+1)/2*gamma * log(M^2/f(M,gamma))
    end

end

