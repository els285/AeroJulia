# Rayleigh flow

"""
Rayleigh flow is frictionless, non-adiabatic flow where heat addition is considered
The stagnation temperature is not constant in this flow
"""

module Rayleigh

    function P_ratio(M,gamma)
        return (gamma+1)/(1+gamma*M^2)
    end

    function Rho_ratio(M,gamma)
        return 1/M^2*1/P_ratio(M,gamma)
    end 

    function T_ratio(M,gamma)
        return P_ratio(M,gamma)^2*M^2
    end

    function V_ratio(M,gamma)
        return P_ratio(M,gamma)*M^2
    end 

    function P_stag_ratio(M,gamma)
        return P_ratio(M,gamma)*(2/(gamma+1)*(1+(gamma-1)/2*M^2))^(gamma/(gamma-1))
    end

    function T_stag_ratio(M,gamma)
        return 2*(gamma+1)*M^2/(1+gamma*M^2)^2*(1+(gamma-1)/2*M^2)
    end

end
