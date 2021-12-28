# Nozzles

module Nozzles

    function AreaRelation(M,gamma)
        return 1/M^2 * (2/(gamma+1) * (1+(gamma-1)/2*M^2))^(gamma+1/gamma-1)
    end

    function AreaRelation_fromFlow(flow)
        return AreaRelation(flow.MachNumber , flow.FluidObj.Ratio_of_Specific_Heats)
    end
    
end 