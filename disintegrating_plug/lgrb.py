import sympy
from probe import show
from equation_of_motion import (alpha, M_1, gamma_i, t_i)
from rhd import c
from rayleigh_taylor import (
    w_0,
    rho_b0,
    derive_fluid_frame_breakup_time,
    derive_breakup_lorentz_factor
    )
import astropy.constants as C
import astropy.units as U

sympy.var('x rho_a omega L M_e R_e', positive=True)

def calc_hydrostatic_density():

    return (M_e/R_e**3)*(x/R_e)**omega

def eruption_criterion():

    return L/(alpha**2*x**2*c) - rho_a*c**2/alpha**4

def calc_eruption_depth():

    _ = eruption_criterion()
    _ = _.subs(rho_a, calc_hydrostatic_density())
    return sympy.solve(_, x)[0]

def calc_eruption_density():

    _ = calc_hydrostatic_density()
    _ = _.subs(x, calc_eruption_depth())
    _ = sympy.expand_power_base(_, force=True)
    _ = _.simplify()
    return _

def calc_eruption_mass():

    _ = calc_eruption_density()
    _ = _*x**3
    _ = _.subs(x, calc_eruption_depth())
    _ = sympy.expand_power_base(_, force=True)
    _ = _.simplify()
    return _

fiducial = {alpha:0.1,
            omega:3,
            L:1e50,
            M_e:10*2e33,
            R_e:3*7e10,
            c:3e10}

def massage(expr):

    _ = expr
    _ = _.subs(w_0, x)
    _ = _.subs(t_i, x/c)
    _ = _.subs(x, calc_eruption_depth())
    _ = _.subs(M_1, calc_eruption_mass())
    _ = _.subs(rho_b0, calc_eruption_density())
    _ = _.subs(gamma_i, 1/alpha)
    _ = _.n()
    return _

def estimate_breakup_radius():

    args = (L,M_e,R_e,c,w_0)

    _ = c*derive_fluid_frame_breakup_time()
    _ = massage(_)
    expr = _.subs(omega, 3)
    _ = _.subs(fiducial)
    f = sympy.lambdify(args, _)
    value = f(*(itm.subs(fiducial) for itm in args))
    return [expr, value, value/7e10*sympy.Symbol(r'R_{\odot}')]

def estimate_breakup_lorentz_factor():

    args = (L,M_e,R_e,c,w_0)

    _ = derive_breakup_lorentz_factor()
    _ = massage(_)
    _ = _.n()
    expr = _.subs(omega, 3)
    _ = _.subs(fiducial)
    f = sympy.lambdify(args, _)
    value = f(*(itm.subs(fiducial) for itm in args))
    return [expr, value]

if __name__ == '__main__':

    show(locals())
