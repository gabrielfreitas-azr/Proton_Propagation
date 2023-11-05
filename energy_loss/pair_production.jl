include("../physical_constants.jl") 
using QuadGK

function beta_pair(E)

    if E > 1E23
        return 0
    end
    
    lorentz_factor =  E / joule_to_ev(m_p * c^2)
    nu = joule_to_ev(m_e * c^2) / (2 * lorentz_factor * k_B * T)

    return pair_const * f_nu(nu) * 1 / E
end

function b0_pair(E)

    return E * beta_pair(E)
end

function db0_pair_dE(E)

    h = 1e-6 * E

    return (b0_pair(E + h) - b0_pair(E - h)) / (2 * h)
end

function phi(xi)

    if xi < 25
        return pi/12 * (xi - 2)^4 / (1 + 0.8048*(xi - 2) + 0.1459*(xi - 2)^2 + 1.137E-3*(xi - 2)^3 - 3.879E-6*(xi - 2)^4)
    else
        return xi*(-86.7 + 50.96*log(xi) - 14.45*log(xi)^2 + 8/3*log(xi)^3)/(1 - (2.910*xi^(-1) + 78.35*xi^(-2) + 1837*xi^(-3)))
    end
end

function f_nu(nu)

    return nu ^ 2 * quadgk(xi -> phi(xi) / (exp(xi * nu) - 1), 2, Inf)[1]
end