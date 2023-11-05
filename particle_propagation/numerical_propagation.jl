include("../energy_loss/energy_loss.jl")
include("einstein_de_sitter.jl")
include("runge_kutta_4.jl")

function f(z, E)

    return -E * beta(E, z) * dt_dz(z)
end

function proton_propa(E_src, z_src)

    return runge_kutta_4(f, z_src, E_src, 0, 1000)
end

function inv_proton_propa(E, z_src)

    return runge_kutta_4(f, 0, E, z_src, 1000)
end

function E_g(E, z_src)

    return inv_proton_propa(E, z_src)[2][end]
end
