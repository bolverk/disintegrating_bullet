import sympy
from probe import show
from rhd import (r, c)
from equation_of_motion import (alpha, M_1, gamma_i, t_i)
from rayleigh_taylor import (
    w_0,
    rho_b0,
    derive_fluid_frame_breakup_time,
    derive_breakup_lorentz_factor
    )

sympy.var('h_e M_e L', positive=True)

def make_density_profile():

    return (M_e/r**3)*sympy.exp(-r**2/h_e**2)

def make_eruption_radius():

    return h_e*sympy.sqrt(sympy.log(M_e*c**3/L/alpha**2/h_e))

def make_eruption_density():

    return alpha**2*L/(h_e**2*c**3)

def make_eruption_mass():

    return alpha**4*L*h_e/c**3

def massage(expr):

    _ = expr
    _ = _.subs(w_0, r)
    _ = _.subs(t_i, r/c)
    _ = _.subs(r, make_eruption_radius())
    _ = _.subs(M_1, make_eruption_mass())
    _ = _.subs(rho_b0, make_eruption_density())
    _ = _.subs(gamma_i, 1/alpha)
    _ = _.n()
    return _

fiducial = {alpha:0.1,
            L:1e50,
            M_e:1e-3*2e33,
            h_e:0.008*7e10,
            c:3e10}

def estimate_breakup_radius():

    _ = c*derive_fluid_frame_breakup_time()
    _ = massage(_)
    expr = _
    _ = _.subs(fiducial)
    value = _
    return [expr, value, value/7e10*sympy.Symbol('R_{\odot}')]

if __name__ == '__main__':

    show(locals())

