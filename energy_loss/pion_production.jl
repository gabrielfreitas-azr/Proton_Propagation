include("../physical_constants.jl")

function b0_pion(E)
    
    return E * beta_pion(E)
end

function db0_pion_dE(E)

    h = 1e-6 * E

    return (b0_pion(E + h) - b0_pion(E - h)) / (2 * h)
end

function beta_pion(E)

    if E < 1e19
        return 0
    end

    E20 = E * 10^(-20)

    return c * 3.2408E-23 / 3.17098E-8 * exp(-4 / E20) / 13.6 * (1 + 4 / E20 + 0.5 * (4 / E20)^2)
end
