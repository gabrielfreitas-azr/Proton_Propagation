include("pair_production.jl")
include("pion_production.jl")
include("adiabatic_loss.jl")

function beta(E, z)
    
    return (1 + z)^3 * beta_pair(E * (1 + z)) + 
           (1 + z)^3 * beta_pion(E * (1 + z)) + 
           (1 + z)^(3/2) * beta_rsh()
end

function b(E, z)

    return (1 + z)^2 * b0(E * (1 + z))
end

function b0(E)

    return b0_pair(E) + b0_pion(E)
end

function db0_dE(E)

    h = 1e-6 * E

    return (b0(E + h) - b0(E - h)) / (2 * h)
end