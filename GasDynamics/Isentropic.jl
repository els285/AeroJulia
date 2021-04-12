module Isentropic

mutable struct IsentropicRatios
    Temperature_Ratio::Float64
    Pressure_Ratio::Float64
    Density_Ratio::Float64
end

function T(M1,M2,gamma)

    return (1 + 0.5*(gamma-1)*M1^2)  /  (1 + 0.5*(gamma-1)*M2^2)

end

function P(M1,M2,gamma)

    return T(M1,M2,gamma)^(gamma/(gamma-1))

end

function Rho(M1,M2,gamma)

    return T(M1,M2,gamma)^(1/(gamma-1))

end


function Generate_IsentropicRatios(M1,M2,gamma)

    Temperature_Ratio   = T(M1,M2,gamma)
    Pressure_Ratio      = P(M1,M2,gamma)
    Density_Ratio       = Rho(M1,M2,gamma)

    return IsentropicRatios(Temperature_Ratio,Pressure_Ratio,Density_Ratio)

end

end