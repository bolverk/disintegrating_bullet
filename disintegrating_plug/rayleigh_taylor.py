import sympy
from probe import show
from caching import memory

from rhd import c, p, eta, t, nu
from equation_of_motion import (
    p_i, t_i, gamma_i, gamma,
    calc_pressure_history,
    calc_asymptotic_gamma)

alpha = sympy.Symbol('alpha')
rho_b0 = sympy.Symbol('rho_b0', positive=True)

@memory.cache
def derive_instability_growth():

    tprime = sympy.Symbol("t'", positive=True)
    aprime = sympy.Symbol("a'", positive=True)
    sympy.var('k rho_b w w_1 rho_b1', positive=True)

    _ = tprime**2*k*aprime*rho_b*c**2/p
    _ = _.subs(rho_b, rho_b1*(p/p_i)**(1/eta))
    _ = _.subs(rho_b1, gamma_i*rho_b0)
    _ = _.subs(k, 1/w)
    _ = _.subs(w, w_1*(p_i/p)**(1/eta)*(t_i/t)**2)
    _ = _.subs(aprime, c*gamma(t)/t)
    _ = _.subs(tprime, t/gamma(t))
    _ = _.subs(gamma(t), calc_asymptotic_gamma())
    _ = _.subs(p, calc_pressure_history())
    return _

@memory.cache
def derive_fluid_frame_breakup_time():

    _ = derive_instability_growth()
    _ = _.subs(p_i, rho_b0*c**2*gamma_i**4)
    _ = _.subs({eta:sympy.Rational(4,3),
                nu:2})
    _ = _.simplify()
    _ = sympy.solve(_ - 1, t)[0]
    _ = sympy.expand_power_base(_, force=True)
    _ = _.simplify()
    return _

def derive_breakup_lorentz_factor():

    _ = calc_asymptotic_gamma()
    _ = _.subs({nu:2, eta:sympy.Rational(4,3)})
    _ = _.subs(t, derive_fluid_frame_breakup_time())
    _ = _.subs(p_i, rho_b0*c**2*gamma_i**4)
    _ = sympy.expand_power_base(_, force=True)
    _ = _.simplify()
    return _

if __name__ == '__main__':

    show(locals())
