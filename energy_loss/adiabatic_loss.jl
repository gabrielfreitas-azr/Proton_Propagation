include("../physical_constants.jl")

function b0_rsh(E)
    
    return E * beta_rsh()
end

function db0_rsh_dE()
    
    return beta_rsh()
end

function beta_rsh()

    return H_0
end
