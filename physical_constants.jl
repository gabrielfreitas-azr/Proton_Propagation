#
#   Definition of useful physical constants.
#

global c = 2.99792458e8         # m/s
global h_bar = 6.6e-16          # eV.s
global m_e = 9.1093837015e-31   # Kg
global m_p = 1.67262192369e-27  # Kg
global k_B = 8.6e-5             # eV/K
global alpha = 0.0072973525693  # -
global e = 1.602176634e-19      # C
global T = 2.73                 # K
global H_0 = 75 *1.022e-12      # 1/yr 
global r_0 = 2.82e-15           # m

function joule_to_ev(joule)

    return joule * 6.242e18
end

#   Pair-production energy rate constant.

global pair_const = (alpha * r_0^2 * (joule_to_ev(m_e * c^2) * k_B * T)^2 * c)/
       (pi^2 * h_bar^3 * c^3) / 3.17098e-8 # eV/yr