import sympy
from probe import show
from caching import memory

from rhd import c, p, eta, t
from equation_of_motion import (
    p_i, gamma,
    calc_pressure_history,
    calc_asymptotic_gamma)

@memory.cache
def derive_instability_growth():

    tprime = sympy.Symbol("t'", positive=True)
    aprime = sympy.Symbol("a'", positive=True)
    sympy.var('k rho_b w w_1 varpi', positive=True)

    _ = tprime**2*k*aprime*rho_b*c**2/p
    _ = _.subs(k, 1/w)
    _ = _.subs(w, w_1*(p_i/p)**(1/eta)*(varpi/c/t)**2)
    _ = _.subs(aprime, c*gamma(t)/t)
    _ = _.subs(tprime, t/gamma(t))
    _ = _.subs(gamma(t), calc_asymptotic_gamma())
    _ = _.subs(p, calc_pressure_history())
    return _

if __name__ == '__main__':

    show(locals())
