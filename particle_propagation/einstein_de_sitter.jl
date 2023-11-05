include("../physical_constants.jl")

function z(t_)

    return (t_0 / t_)^(2/3) - 1
end

function t(z_)
       
    return (2 / 3 * 1 / H_0) * (1 + z_)^(-3/2)
end



function dt_dz(z)

    return -1 / H_0 * (1 + z)^(-5/2)
end

global t_0 = t(0)

function H(z)
    return H_0 * (1 + z)^(3/2)
end

